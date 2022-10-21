
import os
import subprocess
import pandas as pd

import more_itertools as mit

#list of bases recorded at each position that meet a minimum qscore. can then compare against genome later. 1-index rather than 0
def positionalData(*, sortedBamfile, outfile, qvalCutoff):
    tempFile = sortedBamfile.split('.')[0] + '_temp.sam'
    output = {}
    command  = ' '.join(['samtools view -f 0x1 -F 0x900', sortedBamfile,">", tempFile])
    os.system(command)
    with open(tempFile, 'r') as infile:
        while True:
            line1 = infile.readline()
            line2 = infile.readline()
            if not line2:
                break
            firstMate = line1.split('\t')
            secondMate = line2.split('\t')
            if firstMate[0] != secondMate[0]:
                print('Odd number of matepairs in file...fix')
                break
            #for each, store mapped positions (and small deletions, not counting them) in one list, and deletions in another. Then compare.
            #ignore chimeras
            if firstMate[2] == secondMate[2]:
                currentSeq = {}
                segment = firstMate[2]
                for read in [firstMate, secondMate]:
                    CIGAR = read[5]
                    currentPosSeq = int(read[3])
                    score = read[10]
                    sequence = read[9]
                    #start positions and mapping. Just go with lowest score at position rather than handling in a Bayesian framework. There is a lot of discussion about that.
                    CIGARval = []
                    curstring = ''
                    for character in CIGAR:
                        if not character.isalpha():
                            curstring += character
                        else:
                            CIGARval.append((int(curstring), character))
                            curstring = ''
                    currentPosRead = 1
                    #For D, N, advance position in sequence but not read. For I, S, advance position in read but not sequence
                    #For M, add in
                    for value in CIGARval:
                        if value[1] == 'M':
                            #add all mapped positions to a list
                            basesAdded = 0

                            while basesAdded < value[0]:
                                if currentPosSeq not in currentSeq or ord(score[currentPosRead-1]) - 33 > currentSeq[currentPosSeq][1]:
                                    currentSeq[currentPosSeq] = (sequence[currentPosRead-1], ord(score[currentPosRead-1]) - 33)
                                basesAdded += 1
                                currentPosSeq += 1
                                currentPosRead +=1
                        elif value[1] == 'N' or value[1] == 'D':
                            currentPosSeq += value[0]
                        elif value[1] == 'I' or value[1] == 'S':
                            currentPosRead += value[0]
                for position in currentSeq:
                    if segment not in output:
                        output[segment] = {}
                    if position not in output[segment]:
                        output[segment][position] = {'A':0, 'G':0,'C':0,'T':0}
                    if currentSeq[position][1] >= qvalCutoff:
                        #for any reasonable qValCufoff none should be N's, but just to be anal retentive
                        base = currentSeq[position][0]
                        if base != 'N':
                            output[segment][position][base] += 1
        with open(outfile, 'w') as outfile:
            #might as well put it all in order in a tsv
            outfile.write('\t'.join(['segment','position','A','G','C','T']) + '\n')
            for segment in sorted(output.keys()):
                for position in sorted(output[segment]):
                    writeline = '\t'.join([segment, str(position)])
                    for base in ['A','G','C','T']:
                        writeline += '\t' + str(output[segment][position][base])
                    writeline += '\n'
                    outfile.write(writeline)
        arguments = f"rm -f {tempFile}"






def reduceComappedToRetained(rawCounts, mapped, curatedMappedFile):
    rawCounts = pd.read_csv(rawCounts, sep='\t')
    mapped = pd.read_csv(mapped, sep='\t')
    rawCounts['combined'] = rawCounts.fivePrime.astype(str) + ':' + rawCounts.threePrime.astype(str)
    #first pass
    mapped['combined1'] = mapped.fivePrime_1.astype(str) + ':' +  mapped.threePrime_1.astype(str)
    mapped['combined2'] = mapped.fivePrime_2.astype(str) + ':' +  mapped.threePrime_2.astype(str)

    mapped = mapped.loc[mapped.combined1.isin(rawCounts.combined)]
    mapped = mapped.loc[mapped.combined2.isin(rawCounts.combined)]
    mapped = mapped.drop(columns=['combined1', 'combined2'])
    mapped.to_csv(curatedMappedFile, sep='\t', index=False)



'''raw counts for fragments that overlap a given deletion junction. This does NOT consider segments with more than 1 deletion.
This script does produce an extra logfile recording their existance....to be handled
downstream by an additional script. We can NEVER technically rule out all extra deletions, but do want to handle those that
are observed. Provide a name-sorted bamfile. Output chimeras as well '''
def rawCounts(*, junctionFile, sortedBamfile, outDirectory, requiredMapped):
    #load in Junctions
    junctions = pd.read_csv(junctionFile, sep='\t', names= ['segment', 'fivePrime','threePrime','polarity'])
    junctions.drop(columns=['polarity'], inplace=True)
    junctions['counts'] = 0
    multiDel = pd.DataFrame(columns=['segment', 'fivePrime_1', 'threePrime_1', 'fivePrime_2', 'threePrime_2'])
    #STAR uses intron positions, I prefer keeping the last mapped position, so adjust here
    junctions['fivePrime'] = junctions.fivePrime.astype(int) - 1
    junctions['threePrime'] = junctions.threePrime.astype(int) + 1
    sampleName = sortedBamfile.split('/')[-1].split('.')[0]
    junctionCounts = '/'.join([outDirectory,sampleName + '_fragmentCountsRaw.tsv'])
    ChimeraCounts = '/'.join([outDirectory,sampleName + '_chimeraCountsRaw.tsv'])
    multiDeletions = '/'.join([outDirectory, sampleName + '_coMappedDeletions.tsv'])
    #if not (os.path.isfile(junctionCounts) and os.path.isfile(multiDeletions)):
    if not (os.path.isfile(junctionCounts) and os.path.isfile(multiDeletions)):
        #sort so we can make sure we do this as per fragment rather than per read. Then pull ONLY concordantly mapped fragments wherein
        #both reads mapped. Lastly, don't waste time on uninformative reads, for the most part we are only going to look at shifts within deletions.
        #One caveat, we cull ALL multimapped reads rather than assign them using some Bayesian framework. Maybe worth doing later, maybe not. If absolute
        #numbers are required, a Bayesian approach would be appropriate, but the relative numbers here seems ok to use as-is.
        #Also, only pull informative fragments (ie those wherein at least one read maps to a deletion junction)
        tempFile = sortedBamfile.split('.')[0] + '_temp.sam'
        command  = ' '.join(['samtools view -f 0x1 -F 0x900', sortedBamfile, "| awk \'{line=$0; m=$6; getline; if((m ~ /N/) || ($6 ~ /N/)) {print line; print;}}\' > ", tempFile])
        os.system(command)
        with open(tempFile, 'r') as infile:
            while True:
                line1 = infile.readline()
                line2 = infile.readline()
                if not line2:
                    break
                firstMate = line1.split('\t')
                secondMate = line2.split('\t')
                if firstMate[0] != secondMate[0]:
                    print('Odd number of matepairs in file...fix')
                    break
                #for each, store mapped positions (and small deletions, not counting them) in one list, and deletions in another. Then compare.
                #ignore chimeras
                if firstMate[2] == secondMate[2]:
                    segment = firstMate[2]
                    #if neither read has a deletion junction ignore
                    deletions = []
                    mapped = []
                    #start positions and mapping
                    for CIGAR in [(firstMate[5], int(firstMate[3])), (secondMate[5], int(secondMate[3]))]:
                        CIGARval = []
                        curstring = ''
                        for character in CIGAR[0]:
                            if not character.isalpha():
                                curstring += character
                            else:
                                CIGARval.append((int(curstring), character))
                                curstring = ''
                        currentPos = int(CIGAR[1])
                        for value in CIGARval:
                            if value[1] == 'M':
                                #add all mapped positions to a list
                                mapped += list(range(currentPos, currentPos + value[0]))
                                currentPos += value[0]
                            elif value[1] == 'N':
                                deletions += list(range(currentPos, currentPos + value[0]))
                                currentPos += value[0]
                            else:
                                currentPos += value[0]
                    mapped = set(mapped)
                    deletions = set(deletions)
                    if len(deletions) != 0:
                        #no conflict regarding deletions and mapped...or toss
                        if len(mapped & deletions) == 0:
                            #now just see if, between the two reads, we get at least requiredMapped mapped across any given junction.
                            deletionJunctions = []
                            deletions = [list(group) for group in mit.consecutive_groups(sorted(list(deletions)))]
                            for deletion in deletions:
                                #first and last base in a group
                                fivePrime = deletion[0]
                                threePrime = deletion[-1]
                                #see if at least number of mapped reads on each side
                                if len(set(range(fivePrime-requiredMapped, fivePrime)) & set(mapped)) == requiredMapped:
                                    if len(set(range(threePrime+1, threePrime+1 + requiredMapped)) & set(mapped)) == requiredMapped:
                                        #junctions are set as last included base, not first excluded
                                        deletionJunctions += [(fivePrime-1, threePrime+1)]
                                        junctions.loc[(junctions.segment == segment) & (junctions.fivePrime == fivePrime-1) & (junctions.threePrime == threePrime+1), 'counts'] += 1
                                deletionJunctions = sorted(deletionJunctions)
                            if len(deletionJunctions) > 1:
                                for deletion in deletionJunctions:
                                    while len(deletionJunctions) > 0:
                                        compare = deletionJunctions.pop(0)
                                        for deletion in deletionJunctions:
                                            multiDel = multiDel.append(pd.DataFrame({'fivePrime_1':[compare[0]], 'threePrime_1':[compare[1]], 'fivePrime_2':[deletion[0]], 'threePrime_2':[deletion[1]] , 'segment':segment}))


        arguments = f"rm -f {tempFile}"
        subprocess.run(arguments, shell=True, check=True)
        junctions.to_csv(junctionCounts, sep='\t', index=False)
        #collapse the multimaps to unique instances, but go ahead and include counts
        multiDel['counts'] = 1
        multiDel = multiDel.groupby(['segment', 'fivePrime_1', 'threePrime_1', 'fivePrime_2', 'threePrime_2']).counts.sum().reset_index()
        multiDel.to_csv(multiDeletions, sep='\t', index=False)
    reduceComappedToRetained(junctionCounts, multiDeletions, multiDeletions.split('.')[0] + 'Curated.tsv')







'''Go back through all informative reads, looking only at comapped reads. Record all deletions, including triple. To do so,
consider all possible pairwise and 3-wise combinations when calculating exclusion. '''
def relativeComapping(*,sortedBamfile, junctionFile, outfile, requiredMapped):
    junctions = pd.read_csv(junctionFile, sep='\t')
    junctions['deletion_1_only_counts'] = 0
    junctions['deletion_2_only_counts'] = 0
    junctions['deletion_counts'] = 0
    if not (os.path.isfile(outfile)):
        #sort so we can make sure we do this as per fragment rather than per read. Then pull ONLY concordantly mapped fragments wherein
        #both reads mapped. Lastly, don't waste time on uninformative reads, for the most part we are only going to look at shifts within deletions.
        #One caveat, we cull ALL multimapped reads rather than assign them using some Bayesian framework. Maybe worth doing later, maybe not. If absolute
        #numbers are required, a Bayesian approach would be appropriate, but the relative numbers here seems ok to use as-is.
        #Also, only pull informative fragments (ie those wherein at least one read maps to a deletion junction)
        tempFile = sortedBamfile.split('.')[0] + '_temp.sam'
        command  = ' '.join(['samtools view -f 0x1 -F 0x900', sortedBamfile, "| awk \'{line=$0; m=$6; getline; if((m ~ /N/) || ($6 ~ /N/)) {print line; print;}}\' > ", tempFile])
        os.system(command)
        with open(tempFile, 'r') as infile:
            while True:
                line1 = infile.readline()
                line2 = infile.readline()
                if not line1:
                    break
                if not line2:
                    break
                firstMate = line1.split('\t')
                secondMate = line2.split('\t')
                if firstMate[0] != secondMate[0]:
                    print('Odd number of matepairs in file...fix')
                    break
                #for each, store mapped positions (and small deletions, not counting them) in one list, and deletions in another. Then compare.
                if firstMate[2] == secondMate[2]:
                    segment = firstMate[2]
                    #if neither read has a deletion junction ignore
                    deletions = []
                    mapped = []
                    #start positions and mapping
                    for CIGAR in [(firstMate[5], int(firstMate[3])), (secondMate[5], int(secondMate[3]))]:
                        CIGARval = []
                        curstring = ''
                        for character in CIGAR[0]:
                            if not character.isalpha():
                                curstring += character
                            else:
                                CIGARval.append((int(curstring), character))
                                curstring = ''
                        currentPos = int(CIGAR[1])
                        for value in CIGARval:
                            if value[1] == 'M':
                                #add all mapped positions to a list
                                mapped += list(range(currentPos, currentPos + value[0]))
                                currentPos += value[0]
                            elif value[1] == 'N':
                                deletions += list(range(currentPos, currentPos + value[0]))
                                currentPos += value[0]
                            else:
                                currentPos += value[0]
                    mapped = set(mapped)
                    deletions = set(deletions)
                    if len(deletions) != 0:
                        #no conflict regarding deletions and mapped...or toss
                        if len(mapped & deletions) == 0:
                            #now just see if, between the two reads, we get at least requiredMapped mapped across any given junction.
                            deletionJunctions = []
                            deletionSet = [list(group) for group in mit.consecutive_groups(sorted(list(deletions)))]
                            for deletion in deletionSet:
                                #first and last base in a group
                                fivePrime = deletion[0]
                                threePrime = deletion[-1]
                                #see if at least number of mapped reads on each side
                                if len(set(range(fivePrime-requiredMapped, fivePrime)) & set(mapped)) == requiredMapped:
                                    if len(set(range(threePrime+1, threePrime+1 + requiredMapped)) & set(mapped)) == requiredMapped:
                                        #junctions are set as last included base, not first excluded
                                        deletionJunctions += [(fivePrime-1, threePrime+1)]
                                deletionJunctions = sorted(deletionJunctions)
                            if len(deletionJunctions) >= 1:
                                for deletion in deletionJunctions:
                                    #check everything to see if consistent with known deletions. Deal with appropriately. Then add definitively mapped
                                    for index, row in junctions.loc[(junctions.fivePrime_1 ==deletion[0]) & (junctions.threePrime_1 == deletion[1]
                                            ) & (junctions.segment == segment)].iterrows():
                                        specDel = set(range(row.fivePrime_2 + 1, row.threePrime_2)) | set(range(deletion[0] + 1, deletion[1]))
                                        specMapped = (set(range(row.fivePrime_2 - requiredMapped, row.fivePrime_2 + 1)) | set(range(row.threePrime_2, row.threePrime_2 + requiredMapped)))
                                        if ((len(mapped & specDel) >= 1) | (len(deletions & specMapped) >= 1)):
                                            junctions.loc[index, 'deletion_1_only_counts'] += 1
                                    #now if it was deletion_2
                                    for index, row in junctions.loc[(junctions.fivePrime_2 ==deletion[0]) & (junctions.threePrime_2 == deletion[1]
                                            )  & (junctions.segment == segment)].iterrows():
                                        specDel = set(range(row.fivePrime_1 + 1, row.threePrime_1)) | set(range(deletion[0] + 1, deletion[1]))
                                        specMapped = (set(range(row.fivePrime_1 - requiredMapped, row.fivePrime_1 + 1)) | set(range(row.threePrime_1, row.threePrime_1 + requiredMapped)))
                                        if ((len(mapped & specDel) >= 1) | (len(deletions & specMapped) >= 1)):
                                            junctions.loc[index, 'deletion_2_only_counts'] += 1
                                if len(deletionJunctions) >= 2:
                                    while len(deletionJunctions) > 0:
                                        compare = deletionJunctions.pop(0)
                                        for deletion in deletionJunctions:
                                            junctions.loc[(junctions.fivePrime_1 == compare[0]
                                             ) & (junctions.threePrime_1 == compare[1]) & (junctions.fivePrime_2 == deletion[0]) & (junctions.threePrime_2 == deletion[1]
                                             ) & (junctions.segment == segment), 'deletion_counts'] += 1


        arguments = f"rm -f {tempFile}"
        subprocess.run(arguments, shell=True, check=True)
        junctions.to_csv(outfile, sep='\t', index=False)



'''Originally tried to assign reads by probability, but with sparse data this feels...not great. For now just collapse everything
to most abundant-ie >50% dual junction, single, or triple'''
def finalJunctions(*, rawJunctions, coMapped, outfile, fasta):
    lengths = []
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
                lengths += [pd.DataFrame({'segment':[name], 'length':[length]})]
        lengths = pd.concat(lengths)
    deletionInfo = {}
    coMappedFrame = pd.read_csv(coMapped, sep='\t')
    #just store a list of fraction junctions mapping, will collapse to max for assignment
    for index, row in coMappedFrame.iterrows():
        if row.deletion_counts != 0:
            deletion1 = ':'.join([row.segment, str(row.fivePrime_1), str(row.threePrime_1)])
            deletion2 = ':'.join([row.segment, str(row.fivePrime_2), str(row.threePrime_2)])
            probability1 = float(row.deletion_1_only_counts)/(row.deletion_1_only_counts + row.deletion_counts)
            probability2 = float(row.deletion_2_only_counts)/(row.deletion_2_only_counts + row.deletion_counts)
            if deletion1 in deletionInfo.keys():
                deletionInfo[deletion1] += [probability1]
            else:
                deletionInfo[deletion1] = [probability1]
            if deletion2 in deletionInfo.keys():
                deletionInfo[deletion2] += [probability2]
            else:
                deletionInfo[deletion2] = [probability2]
    multiProb = []
    for deletion in deletionInfo:
        prob = max(deletionInfo[deletion])
        multiProb += [pd.DataFrame({'segment':[deletion.split(':')[0]], 'deletion_junction':[':'.join(deletion.split(':')[1:])],
                        'multiple_probability_max':[prob]})]
    multiProb = pd.concat(multiProb)
    rawCounts = pd.read_csv(rawJunctions, sep='\t')
    rawCounts['deletion_junction'] = rawCounts.fivePrime.astype(str) + ':' + rawCounts.threePrime.astype(str)
    rawCounts = pd.merge(rawCounts, lengths, on='segment')
    rawCounts['deletion_length'] = rawCounts.threePrime - rawCounts.fivePrime - 1
    rawCounts['final_length'] = rawCounts.length  - rawCounts.deletion_length
    rawCounts = pd.merge(rawCounts, multiProb, on=['segment', 'deletion_junction'], how='left').fillna(0)
    rawCounts = rawCounts.drop(columns = ['length'])
    rawCounts.to_csv(outfile, sep='\t', index=False)


'''fragments that align across center. Take only readpairs that do not have any deletion junctions wherein at least one read overlaps the central portion of the gene.
go ahead and make a bamfile consisting of unspliced reads that map the central portion'''
def centralDepth(fasta, sortedBamfile, outfile, maxFragmentSize):
    numbers = {}
    #try and do the longest steps just once rather than for each flu segment, although the argument could be made that splitting would be faster
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
    for segment in lengths:
        center = str(int(lengths[segment]/2))

        command  = ''.join(['samtools view -f 0x1 -F 0x900 ', sortedBamfile,
        " | awk \'{line=$0; m=$6; s=$3; p=$4; l=$9; getline; if((m !~ /N/) && ($6 !~ /N/) && ((($4 + $9 > ",
        center,
        ") && ($4 < ",
        center,
        ")) || ((p + l > ",
        center,
        ") && (p < ",
        center,
        "))) && ((l < ",
        str(maxFragmentSize),
        ") && ($9 < ",
        str(maxFragmentSize),
        ")) && ((s == \"", segment, "\") && ($3 == \"", segment, "\"))) {print line; print;}}\' | wc -l"])
        process = subprocess.Popen(command, shell=True,  stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        #using float as should always be an int, this is just a little check for me
        number = (float(stdout.decode('utf-8')))/float(2)
        numbers[segment] = number
    pd.DataFrame.from_dict(numbers, orient='index', columns=['count']).reset_index().rename(columns={'index':'segment'}).to_csv(outfile, sep='\t', index=False)
