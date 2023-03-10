'''take a blast file produced from bash script delcaller and identifies deletion junctions'''

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import pandas as pd
#Identify blast reads that map across a junction. Should be two (and no more) blast hits across a junction. For repetitive elements (ie, reads map twice at junction), go ahead and take the
#bases prior to the deletion (ie 5'), as it is an academic decision anyway. Nevertheless, want to make sure I only identify one deletion and not two that result in identical sequence

def delMap(blastin, outfile, fasta):
    sequenceDict = {}
    segName = ''
    with open(fasta, 'r') as infile:
        for line in infile:
            if line[0] == '>':
                segName = line[1:-1]
            elif segName != '':
                sequenceDict[segName] = line[:-1]
                segName = ''
    deletionFrame = {}
    with open(blastin) as blastfile:
        line1 = blastfile.readline()
        while(line1 != ""):
            line2 = blastfile.readline()
            lineseg1 = line1.split("\t")
            lineseg2 = line2.split("\t")

            #if across junction, adjacent lines will match one another
            if lineseg1[0] == lineseg2[0]:
                #toss PCR hybrids
                if lineseg1[7] == lineseg2[7]:
                    line3 = blastfile.readline()
                    lineseg3 = line3.split("\t")
                    #If three reads in a row map, rather than deal with ambiguity, just move on. Make sure you capture all the reads that behave that way though, go until a new name is hit
                    if lineseg3[0] == lineseg2[0]:
                        line1 = blastfile.readline()
                        while line1.split('\t')[0] == lineseg2[0]:
                            line1 = blastfile.readline()
                    else:
                        #Definition of blastfile, order is read id, read length, query start, query end, sequence start, sequence end, e-value, segment ID
                        #If all bases unaccounted for (ie not both sides of junction mapping) toss this thing
                        #Lets determine what order the sequences should go in, and their polarity
                        #Order of sequences according to actual read
                        if(int(lineseg1[2]) < int(lineseg2[2])):
                            first = lineseg1
                            second = lineseg2
                        else:
                            first = lineseg2
                            second = lineseg1
                        #Check that the reads actually span across a junction. If there are "missing" bases throw it away
                        if (int(first[3]) >= (int(second[2])-1)):
                            #Ok, polarity check
                            forward = (int(first[4]) < int(first[5]))
                            #Now both reads should have identical polary, or else there is an inversion!
                            if forward == (int(second[4]) < int(second[5])):
                                #So now I know the read order, the polarity, that they map across a junction, and that they do not represent an inversion.
                                #Still need to check for other, bizarre PCR hybrids, but will do so within following area with a toggle
                                writeit = False
                                if forward:
                                    if(int(first[5]) < int(second[4])):
                                        #need to know if there is overlap
                                        overlap = int(first[3]) - int(second[2]) + 1
                                        #need to know how to handle overlap. Which bases are given to which sequence
                                        given_left = 0
                                        given_right = 0
                                        hit_mismatch = False
                                        if overlap > 0:
                                            overlapsequence = first[8][-overlap -1:-1]
                                            delCount = 0
                                            newoverlapsequence = ""
                                            for character in overlapsequence:
                                                if character != '-':
                                                    newoverlapsequence += character
                                                else:
                                                    delCount += 1
                                            overlapsequence = first[8][-overlap -1 - delCount:-overlap - 1] + newoverlapsequence
                                            refSeqLeft = sequenceDict[first[7]][int(first[5])- overlap:int(first[5])]
                                            refSeqRight =  sequenceDict[second[7]][int(second[4]) - 1: int(second[4]) + overlap - 1]
                                            #now look for any mismatches that crept in. As convention, I will always give bases to the 5' end first
                                            #if mismatches on both, base will still be given to left end. Once a mismatch is hit that is present in the
                                            #3' end, all remaining characters are given to that position
                                            for position, character in enumerate(overlapsequence):
                                                if hit_mismatch:
                                                    given_right += 1
                                                elif character == refSeqLeft[position]:
                                                    given_left += 1
                                                elif character == refSeqRight[position]:
                                                    given_right += 1
                                                    hit_mismatch = True
                                                else:
                                                    given_left += 1

                                        junction_left = str(int(first[5]) - overlap + given_left)
                                        junction_right = str(int(second[4]) + overlap - given_right)
                                        writeit = True

                                else:
                                    #Now for opposite polarity
                                    if(int(first[5]) > int(second[4])):
                                        overlap = int(first[3]) - int(second[2]) + 1
                                        given_left = 0
                                        given_right = 0
                                        hit_mismatch = False
                                        refSeqLeft = ""
                                        refSeqRight = ""
                                        if overlap > 0:
                                            overlapsequence = first[8][-overlap -1:-1]
                                            delCount = 0
                                            newoverlapsequence = ""
                                            for character in overlapsequence:
                                                if character != '-':
                                                    newoverlapsequence += character
                                                else:
                                                    delCount += 1
                                            overlapsequence = reverseComplement(first[8][-overlap -1 - delCount:-overlap - 1] + newoverlapsequence)
                                            refSeqLeft = sequenceDict[second[7]][int(second[4]) - overlap:int(second[4])]
                                            refSeqRight = sequenceDict[first[7]][int(first[5]) - 1: int(first[5]) + overlap - 1]
                                            #now look for any mismatches that crept in. As convention, I will always give bases to the 5' end first
                                            #if mismatches on both, base will still be given to left end. Once a mismatch is hit that is present in the
                                            #3' end, all remaining characters are given to that position
                                            for (position, character) in enumerate(overlapsequence):
                                                if hit_mismatch:
                                                    given_right += 1
                                                elif character == refSeqLeft[position]:
                                                    given_left += 1
                                                elif character == refSeqRight[position]:
                                                    given_right += 1
                                                    hit_mismatch = True
                                                else:
                                                    given_left += 1
                                        junction_left = str(int(second[4]) - overlap + given_left)
                                        junction_right = str(int(first[5]) + overlap - given_right)
                                        writeit = True


                                if writeit:
                                    segment = first[7]
                                    if segment in deletionFrame:
                                        if junction_left in deletionFrame[segment]:
                                            if junction_right in deletionFrame[segment][junction_left]:
                                                deletionFrame[segment][junction_left][junction_right] += 1
                                            else:
                                                deletionFrame[segment][junction_left][junction_right] = 1
                                        else:
                                            deletionFrame[segment][junction_left] = {junction_right:1}
                                    else:
                                        deletionFrame[segment] = {junction_left:{junction_right:1}}

                        #Don't accidentally lose line3! Readline advances file!
                        line1 = line3
                else:
                    line1=line2

            else:
                #Don't want to accidentally advance two lines, this ensures you only advance one
                line1 = line2
    with open(outfile, 'w') as outfile:
        outfile.write('\t'.join(['segment', 'left', 'right', 'count']) + '\n')
        for segment in deletionFrame:
            for left in deletionFrame[segment]:
                for right in deletionFrame[segment][left]:
                    outfile.write('\t'.join([segment,left,right, str(deletionFrame[segment][left][right])]) + '\n')




def reverseComplement(sequence):
    sequence = sequence[::-1]
    RC = ""
    for character in sequence:
        if character == 'A':
            RC += 'T'
        elif character == 'G':
            RC += 'C'
        elif character == 'T':
            RC += 'A'
        elif character == 'C':
            RC += 'G'
    return RC









#take the indicated directory containing deletion junction collapsed files, identify deletions that match some threshold (5' OR 3' end)
#compress deletions in matched samples to generate a single list, and append this list to a GTFfile, generating a novel GTF file
def delCombine(GTFfile, samples, outfile, threshold, segments, lengths):
    with open (GTFfile) as GTFin, open(outfile, 'w') as GTFout:
        #deletions will be a list of segments with dictionaries indexed by 5' ends containing 3' ends
        deletions = {}
        for segment in segments:
            deletions[segment] = {}
        for segment in segments:
            for sample in samples:
                for segment in segments:
                    with open(sample + "/" + segment + ".junccollapse") as deletionfile:
                        #throw away header
                        deletionfile.readline()
                        for line in deletionfile:
                            lineseg = re.split('\t|\n', line)
                            if(float(lineseg[3]) >= threshold or float(lineseg[4]) >= threshold):
                                if lineseg[0] in deletions[segment]:
                                    if lineseg[1] not in deletions[segment][lineseg[0]]:
                                        deletions[segment][lineseg[0]].append(lineseg[1])
                                else:
                                    deletions[segment][lineseg[0]] = [lineseg[1]]
        for line in GTFin:
            GTFout.write(line)
        for (length, segment) in enumerate(segments):
            delCount = 1
            for fiveprime in deletions[segment]:
                for threeprime in deletions[segment][fiveprime]:
                    name = segment + "_deletion_" + str(delCount)
                    line1 = "\t".join([segment, "ABR", "exon", '1', fiveprime, ".", "+", ".", 'gene_id "' + name + '"; gene_name "' + name + '"; transcript_id "' + name + '"; tss_id "' + name + '";\n'])
                    line2 = "\t".join([segment, "ABR", "exon", threeprime, lengths[length], ".", "+", ".", 'gene_id "' + name + '"; gene_name "' + name + '"; transcript_id "' + name + '"; tss_id "' + name + '";\n'])
                    GTFout.write(line1)
                    GTFout.write(line2)
                    delCount += 1

#Just the first part of delCombine, returns the appropriate dictionary for other methods
def delCount(GTFfile, samples, outfile, threshold, segments, lengths):
    with open (GTFfile) as GTFin, open(outfile, 'w') as GTFout:
        #deletions will be a list of segments with dictionaries indexed by 5' ends containing 3' ends
        deletions = {}
        for segment in segments:
            deletions[segment] = {}
        for segment in segments:
            for sample in samples:
                for segment in segments:
                    with open(sample + "/" + segment + ".junccollapse") as deletionfile:
                        #throw away header
                        deletionfile.readline()
                        for line in deletionfile:
                            lineseg = re.split('\t|\n', line)
                            if(float(lineseg[3]) >= threshold or float(lineseg[4]) >= threshold):
                                if lineseg[0] in deletions[segment]:
                                    if lineseg[1] not in deletions[segment][lineseg[0]]:
                                        deletions[segment][lineseg[0]].append(lineseg[1])
                                else:
                                    deletions[segment][lineseg[0]] = [lineseg[1]]

#provide dataframe of 1-indexed junctions and convert to 0-indexed exon and intron files. Maintain an "unspliced" flu vRNA for basal mapping using tmo option
#provide list of 1-indexed junctions and convert to 1-indexed intron files. To make first and last base of intron modify slightly.
def junctionsToIntrons(junctionList, fasta, outBase, basesReq):
    lengths = {}
    with open(fasta, 'r') as infile:
        lines = infile.readlines()
        while len(lines) > 0:
            line1 = lines.pop(0)
            line2 = lines.pop(0)
            if line1[0] != '>':
                raise Exception('Incorrectly formatted fasta file, please fix')
            else:
                #don't forget newline characters
                name = line1[1:-1]
                length = len(line2[:-1])
                lengths[name] = length
    with open(outBase + '_exon.tsv', 'w') as exonFile, open(outBase + '_intron.tsv', 'w') as intronFile:
        for segment in lengths:
            exonFile.write('\t'.join([segment, '0', str(lengths[segment]), '+']) + '\n')
        for deletion in junctionList:
            segment, fiveprime, threeprime = deletion.split(':')
            fiveprime = str(int(fiveprime) + 1)
            threeprime = str(int(threeprime) -1)
            if int(fiveprime) >= basesReq - 1:
                if lengths[segment] - int(threeprime) >= basesReq:
                    intronFile.write('\t'.join([segment, fiveprime, threeprime, '+']) + '\n')




def delAnnotate(blastin, outfile, fasta):
    sequenceDict = {}
    segName = ''
    with open(fasta, 'r') as infile:
        for line in infile:
            if line[0] == '>':
                segName = line[1:-1]
            elif segName != '':
                sequenceDict[segName] = line[:-1]
                segName = ''
    with open(blastin) as blastfile, open(outfile, 'w') as outfile:
        line1 = blastfile.readline()
        while(line1 != ""):
            line2 = blastfile.readline()
            lineseg1 = line1.split("\t")
            lineseg2 = line2.split("\t")

            #if across junction, adjacent lines will match one another
            if lineseg1[0] == lineseg2[0]:
                #toss PCR hybrids
                if lineseg1[7] == lineseg2[7]:
                    line3 = blastfile.readline()
                    lineseg3 = line3.split("\t")
                    #If three reads in a row map, rather than deal with ambiguity, just move on. Make sure you capture all the reads that behave that way though, go until a new name is hit
                    if lineseg3[0] == lineseg2[0]:
                        line1 = blastfile.readline()
                        while line1.split('\t')[0] == lineseg2[0]:
                            line1 = blastfile.readline()
                    else:
                        #Definition of blastfile, order is read id, read length, query start, query end, sequence start, sequence end, e-value, segment ID
                        #If all bases unaccounted for (ie not both sides of junction mapping) toss this thing
                        #Lets determine what order the sequences should go in, and their polarity
                        #Order of sequences according to actual read
                        if(int(lineseg1[2]) < int(lineseg2[2])):
                            first = lineseg1
                            second = lineseg2
                        else:
                            first = lineseg2
                            second = lineseg1
                        #Check that the reads actually span across a junction. If there are "missing" bases throw it away
                        if (int(first[3]) >= (int(second[2])-1)):
                            #Ok, polarity check
                            forward = (int(first[4]) < int(first[5]))
                            #Now both reads should have identical polary, or else there is an inversion!
                            if forward == (int(second[4]) < int(second[5])):
                                #So now I know the read order, the polarity, that they map across a junction, and that they do not represent an inversion.
                                #Still need to check for other, bizarre PCR hybrids, but will do so within following area with a toggle
                                writeit = False
                                if forward:
                                    if(int(first[5]) < int(second[4])):
                                        #need to know if there is overlap
                                        overlap = int(first[3]) - int(second[2]) + 1
                                        #need to know how to handle overlap. Which bases are given to which sequence
                                        given_left = 0
                                        given_right = 0
                                        hit_mismatch = False
                                        if overlap > 0:
                                            overlapsequence = first[8][-overlap -1:-1]
                                            delCount = 0
                                            newoverlapsequence = ""
                                            for character in overlapsequence:
                                                if character != '-':
                                                    newoverlapsequence += character
                                                else:
                                                    delCount += 1
                                            overlapsequence = first[8][-overlap -1 - delCount:-overlap - 1] + newoverlapsequence
                                            refSeqLeft = sequenceDict[first[7]][int(first[5])- overlap:int(first[5])]
                                            refSeqRight =  sequenceDict[second[7]][int(second[4]) - 1: int(second[4]) + overlap - 1]
                                            #now look for any mismatches that crept in. As convention, I will always give bases to the 5' end first
                                            #if mismatches on both, base will still be given to left end. Once a mismatch is hit that is present in the
                                            #3' end, all remaining characters are given to that position
                                            for position, character in enumerate(overlapsequence):
                                                if hit_mismatch:
                                                    given_right += 1
                                                elif character == refSeqLeft[position]:
                                                    given_left += 1
                                                elif character == refSeqRight[position]:
                                                    given_right += 1
                                                    hit_mismatch = True
                                                else:
                                                    given_left += 1

                                        junction_left = str(int(first[5]) - overlap + given_left)
                                        junction_right = str(int(second[4]) + overlap - given_right)
                                        writeit = True

                                else:
                                    #Now for opposite polarity
                                    if(int(first[5]) > int(second[4])):
                                        overlap = int(first[3]) - int(second[2]) + 1
                                        given_left = 0
                                        given_right = 0
                                        hit_mismatch = False
                                        refSeqLeft = ""
                                        refSeqRight = ""
                                        if overlap > 0:
                                            overlapsequence = first[8][-overlap -1:-1]
                                            delCount = 0
                                            newoverlapsequence = ""
                                            for character in overlapsequence:
                                                if character != '-':
                                                    newoverlapsequence += character
                                                else:
                                                    delCount += 1
                                            overlapsequence = reverseComplement(first[8][-overlap -1 - delCount:-overlap - 1] + newoverlapsequence)
                                            refSeqLeft = sequenceDict[second[7]][int(second[4]) - overlap:int(second[4])]
                                            refSeqRight = sequenceDict[first[7]][int(first[5]) - 1: int(first[5]) + overlap - 1]
                                            #now look for any mismatches that crept in. As convention, I will always give bases to the 5' end first
                                            #if mismatches on both, base will still be given to left end. Once a mismatch is hit that is present in the
                                            #3' end, all remaining characters are given to that position
                                            for (position, character) in enumerate(overlapsequence):
                                                if hit_mismatch:
                                                    given_right += 1
                                                elif character == refSeqLeft[position]:
                                                    given_left += 1
                                                elif character == refSeqRight[position]:
                                                    given_right += 1
                                                    hit_mismatch = True
                                                else:
                                                    given_left += 1
                                        junction_left = str(int(second[4]) - overlap + given_left)
                                        junction_right = str(int(first[5]) + overlap - given_right)
                                        writeit = True


                                if writeit:
                                    outfile.write('\t'.join([first[0], first[7], junction_left, junction_right]) + '\n')

                        #Don't accidentally lose line3! Readline advances file!
                        line1 = line3
                else:
                    line1=line2

            else:
                #Don't want to accidentally advance two lines, this ensures you only advance one
                line1 = line2

def delAssign(samfile, junctionFile, outFile):
    with open(junctionFile) as junction, open(samfile) as samfile:
        samLine = samfile.readline()
        #every junction read is in the samfile, but not every
        junctionLine = junction.readline()
        output = []
        while junctionLine != '':
            junctionSeg = junctionLine.split('\t')
            samSeg = samLine.split('\t')
            while(samSeg[0] != junctionSeg[0]):
                samSeg = samfile.readline().split('\t')
            #get the cell barcode variable tags so just go through
            #storing everything in memory, so keep that in mind if this ever gets used for a massive dataset
            cellbarcode = ''
            for element in samSeg:
                if element[:5] == 'CB:Z:':
                    cellbarcode = element[5:]
                if element[:5] == 'UB:Z:':
                    umibarcode = element[5:]
            output += [pd.DataFrame({'barcode':[cellbarcode], 'segment':[junctionSeg[1]], 'umibarcode':[umibarcode],
                'fiveprime':[junctionSeg[2]], 'threeprime':[junctionSeg[3]], 'count':[1]})]
            junctionLine = junction.readline()
    output = pd.concat(output)
    #this step might be slow
    #remove same barcode, same UMI, same deletion
    output = output.drop_duplicates().drop(columns=['umibarcode'])
    output = output.groupby(['barcode','segment','fiveprime','threeprime']).sum().reset_index()
    output.to_csv(outFile, sep='\t')

#gives a count of non-split reads at the indicated site of the indicated barcode

def adjReads(deletionFrame, bamfile):
    deletionFrame['nonDeletionCounts'] = 0
    for row in adjReads.iterrows():
        segment = row.segment
        barcode = row.barcode
        threeprime = row.threeprime
        command = 'samtools view ' + bamfile + ' ' + segment + ':' + str(threeprime - 1) ' | awk  awk \'($6 !~ /N/)\' grep -c CB:Z: ' +barcode
        process = subprocess.Popen(command, shell=True,  stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        number = (float(stdout.decode('utf-8')))
        row['nonDeletionCounts'] = number
    return deletionFrame
