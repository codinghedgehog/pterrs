#!/usr/bin/python
#
# Project ERS to ERR lookup

import os
import sys
import argparse
import urllib2
import xml.etree.ElementTree
import re

ENB_URL="https://www.ebi.ac.uk/ena/data/view/{0}&display=xml"
ENB_FILEREPORT_URL="https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession={0}&result=read_run&fields=experiment_accession,run_alias,run_accession"

if __name__ == '__main__':
    argparser = argparse.ArgumentParser()
    argparser.add_argument("project_id", help="Project ID [e.g. PRJEB211])")
    args = argparser.parse_args()

    projectID = args.project_id

    ersDict = {}

    print "Querying {0}".format(ENB_URL.format(projectID))

    # First, retrieve the XML data for the project (primary id) and get the list of ERS ids
    # referenced by it.
    req = urllib2.Request(ENB_URL.format(projectID))
    response = urllib2.urlopen(req)
    projectXML = response.read()
    
    # ERS ids are under /PROJECT_LINKS/PROJECT_LINK/XREF_LINK/<DB = 'ENA_SAMPLE'>/ID
    projectXMLTree = xml.etree.ElementTree.fromstring(projectXML)

    projectERSElements = projectXMLTree.findall(".//PROJECT_LINKS/PROJECT_LINK/XREF_LINK/[DB='ENA-SAMPLE']/ID")

    if len(projectERSElements) == 0:
        print "No ENA-SAMPLE (ERS IDs) found.  Invalid project ID passed?"
        sys.exit(1)
    elif len(projectERSElements) != 1:
        print "Multiple ENA-SAMPLE DB/ID elements found? Not handled in this version."
        sys.exit(1)

    ersString = projectERSElements[0].text

    print "Found these ERS IDs for project {0}: {1}".format(projectID,ersString)

    # ERS data from XML is a comma-delimited value or value range, e.g.
    # ERS005741-ERS005788,ERS011783-ERS012055 
    #
    # Need to expand them out to create the ERS table.

    ersList = ersString.split(',')
    for ersValue in ersList:
        # Single ERS value, just save it off.
        if re.match("^ERS[0-9]{6}$",ersValue):
            ersDict[ersValue] = ""
        # ERS range, expand it out.
        elif re.match("^ERS(?P<startRange>[0-9]{6})-ERS(?P<endRange>[0-9]{6})$", ersValue):
            rangeMatch = re.match("^ERS(?P<startRange>[0-9]{6})-ERS(?P<endRange>[0-9]{6})$", ersValue)
            startRange = rangeMatch.group('startRange')
            endRange = rangeMatch.group('endRange')

            for ersNum in range(int(startRange),int(endRange) + 1):
                ersDict['ERS{:06d}'.format(ersNum)] = ''
        else:
            print "Unrecognized ERS id format: {0} (not singular value nor range?)".format(ersValue)
            sys.exit(1)

    # At this point, have a dictionary of ERS numbers (programmer note: This is what is supposed to 
    # be in the datasheet).
    #
    # Now need to retrieve all the ERX->ERR correlation data from the project
    # The catch is that while ERR is unique, ERX needs to be paired up with the run_alias (which is the
    # file_id in the datasheet) which looks like SC_RUN_5150_1#0 (in datasheet, it is listed as 5150_1_0).
    #
    # But we want to key off of ERX for the lookup, so to make it unique, need a combined key of ERX+RUN_ALIAS
    #
    # NOTE: This may be project-specific, since RUN_ALIAS seems to be free form data.
    #
    req = urllib2.Request(ENB_FILEREPORT_URL.format(projectID))
    response = urllib2.urlopen(req)
    erx_err_text = response.read()

    erxList = []
    erxDict = {}

    # Store this data in a dict (keying off of ERX number and modified run_alias, tying it to ERR number).
    # Data returned is ERX_NUMBER\t\tRUN_ALIAS\tERR_NUMBER, with header, rows delimited by newline
    erx_err_list = erx_err_text.split('\n')
    for erx_err_pair in erx_err_list:
        # Keep it simple (if fragile) -- split on tab \t char to break out ERX, run_alias, and ERR.
        # We'll end up saving off the column headers, too, but that won't matter.
        erx_err_list = erx_err_pair.split('\t')
        # Check length of split, because output may have empty lines
        if len(erx_err_list) == 3:
            ERX_ID = erx_err_list[0]
            RUN_ALIAS = erx_err_list[1]
            ERR_ID = erx_err_list[2]
            erxDict[ERX_ID + RUN_ALIAS] = ERR_ID
            erxList.append(ERX_ID)

    # Remove duplicate ERX numbers.
    erxList = sorted(set(erxList))
    
    # Now can work from ERX to the ERS + RUN_ALIAS (file id) in the datasheet.
    # Can query this information via a lookup of ERX on the EBI ENA site, which will
    # return XML with a list of ERS/File_ID (RUN_ALIAS) pairs, which can then be used
    # in erxDict to find the ERR number associated with it!
    for erxNum in erxList:
        req = urllib2.Request(ENB_URL.format(erxNum))
        response = urllib2.urlopen(req)
        erxXML = response.read()

        # Extract all the related ERS and run_alias (member_name attribute in XML)
        # ERS ids are under /PROJECT_LINKS/PROJECT_LINK/XREF_LINK/<DB = 'ENA_SAMPLE'>/ID
        erxXMLTree = xml.etree.ElementTree.fromstring(erxXML)
        ersElements = erxXMLTree.findall(".//EXPERIMENT/DESIGN/SAMPLE_DESCRIPTOR/POOL/MEMBER")

        # Each MEMBER element should have an accession attribute (ERS) and member_name (RUN_ALIAS/FILE_ID).
        # Use that to create the lookup key to query in the erxDict to find the ERR ID!
        for ersElement in ersElements:
            ersAccession = ersElement.get('accession')
            ersMemberName = ersElement.get('member_name')

            # Now form the look up key - ERX (not ERS!) + RUN_ALIAS - to get the associated ERR.
            # NOTE: This part is project specific, depending on the RUN_ALIAS formatting.
            erxDictKey = erxNum + 'SC_RUN_' + ersMemberName
            try:
                if erxDict[erxDictKey]:
                    # For output, make member_name match the datasheet's File ID (no SC_RUN and 
                    # convert # to _ in the RUN_ALIAS).
                    fileID = re.sub('SC_RUN_','',ersMemberName)
                    fileID = re.sub('#','_',fileID)
                    print "{0}, File ID {1}: {2}".format(ersAccession, fileID, erxDict[erxDictKey])
            except KeyError:
                print "*** WARNING: No ERR found for {0} with file ID {1}".format(ersAccession,fileID)

