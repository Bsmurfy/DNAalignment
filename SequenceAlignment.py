'''
Brendan Murphy

A program to globally align multiple DNA sequences from a FASTA file.
'''

def main():
    match = 5
    mismatch = -4
    gap = -11

    # read the FASTA file
    fastafull = readfasta('MidCS4.fasta')
    numofseq = min(4, len(fastafull)) # allows up to four DNA sequences, for now

    dnalist = []
    for i in range(0, numofseq):
        dnasequence = fastafull[i][2]
        dnalist.append(dnasequence)
    # run alignment loop until only one sequence is left
    while (numofseq > 1):
        candidates = [] # initialize/clear candidates list
        # make list of DNA sequences from FASTA file
        for seq1 in dnalist:
            for seq2 in dnalist:
                if (dnalist.index(seq2) > dnalist.index(seq1)):
                    candidates = align(candidates, seq1, seq2, match, mismatch, gap)
        scores = []
        # form list of scores to determine max (best) alignment score
        for candidate in candidates:
            scores.append(int(candidate[1]))
        bestscore = max(scores)
        print('test')
        for candidate in candidates:
            if (candidate[1] == bestscore):
                dnalist.insert(0, candidate[0]) # add the aligned sequence to the list of DNA...
                dnalist.remove(candidate[2]) # ...and remove the two old strings
                dnalist.remove(candidate[3])
        numofseq = numofseq - 1

    bestalign = ''.join(dnalist)
    print('The best alignment of these sequences is: ')
    print(bestalign)
            
def align(candidatelist, a, b, match, mismatch, gap):
    m = len(a)
    n = len(b)
    matrix = [[0 for x in range(m)] for y in range(n)]
    # create scoring matrix reflecting match/mismatch scores
    # score_matrix = [[match, mismatch, mismatch, mismatch],
                    #[mismatch, match, mismatch, mismatch],
                    #[mismatch, mismatch, match, mismatch],
                    #[mismatch, mismatch, mismatch, match]]

    # initialize matrix with 0's
    for i in range(0, m):
        for j in range (0, n):    
            matrix[i][j] = 0
    # intialize score as 0
    score = 0
    
    # build matrix (except for first row and first column, which stays 0's)
    if isinstance(a, list):
        for item in a: # if the first item is a list (i.e. has already been assigned)
            for i in range(2, m):
                for j in range(2, n):
                    # note: this uses if/elif and comparison to determine match/mismatch
                    # rather than a scoring matrix
                    if (item[i] == b[j]) and (item[i] != '_'):
                        matrix[i][j] = max(matrix[i-1][j-1] + match,
                                            matrix[i-1][j] + gap,
                                            matrix[i][j-1] + gap)
                    else:
                        matrix[i][j] = max(matrix[i-1][j-1] + mismatch,
                                            matrix[i-1][j] + gap,
                                            matrix[i][j-1] + gap)
            score = score + matrix[m-1][n-1]
    else:
        for i in range(1, m):
            for j in range(1, n):
                # note: this uses if/elif and comparison to determine match/mismatch
                # rather than a scoring matrix
                if (a[i] == b[j]) and (a[i] != '_'):
                    matrix[i][j] = max(matrix[i-1][j-1] + match,
                                        matrix[i-1][j] + gap,
                                        matrix[i][j-1] + gap)
                else:
                    matrix[i][j] = max(matrix[i-1][j-1] + mismatch,
                                        matrix[i-1][j] + gap,
                                        matrix[i][j-1] + gap)
        score = matrix[m-1][n-1]
        
    score = str(score) # the score is the last value in the matrix
    
    i = m - 1
    j = n - 1
    newa = ''
    newb = ''
    # backtrack through the matrix we built
    while (i > 1) and (j > 1): # exits when either sequence is completed
        if (((matrix[i][j] - match) == matrix[i-1][j-1]) or
            ((matrix[i][j] - mismatch) == matrix[i-1][j-1])):
            newa = a[i-1].join(newa) # we're going through the sequences backwards, so add to front of string
            newb = b[j-1].join(newb)
            i = i - 1 # increment counters
            j = j - 1
        elif ((matrix[i][j] - gap) == matrix[i][j-1]):
            newa = '_'.join(newa) # add a gap to one of the strings
            newb = b[j-1].join(newb) # add nucleotide to the other string
            j = j - 1 # only increment one of the counters in gaps
        elif ((matrix[i][j] - gap) == matrix[i-1][j]):
            newa = a[i-1].join(newa)
            newb = '_'.join(newb)
            i = i - 1
    # now, fill the other sequence with gaps as needed
    if (i > 1): 
        while (i > 1):
            newa = a[i-1].join(newa)
            newb = '_'.join(newb)
            i = i - 1
    elif (j > 1):
        while (j > 1):
            newa = '_'.join(newa)
            newb = b[j-1].join(newb)
            j = j - 1
            
    alignedinfo = []
    newseqlist = [newa, newb]
    alignedinfo.extend((newseqlist, score, a, b))
    candidatelist.append(alignedinfo)
    return candidatelist

'''
parseHeader - split out the label from the header line
Parameter: a string starting with ">" and ending without a newline
Return: 
  1. the first item in the string, after the ">", up to the first space
'''
def parseHeaderLine(line):
    header = line[1:]

    label = header.split(' ')[0]

    return label

'''
readfasta - the subroutine that reads the fasta file
Parameter: a filename that must be in fasta format.  The file is
assumed to have:
1. the first line a header line
2. arbitrary blank lines
3. every line (especially including the last) is terminated by a newline
   terminator (enter key)
4. no line has only spaces on it
Return: a list of lists. Each inner list will have three elements:
1. the sequence identifier, the characters between the leading ">"
   and the first space
2. the entire header, the entire first line not including the ">"
3. the sequence, a single string of all the letters with no line terminators
'''
def readfasta(filename):
    resultList = []
    infile = open(filename, 'r')

    # process the first line, which must be a header line
    line = infile.readline()
    headerLine = line.rstrip()

    label = parseHeaderLine(headerLine)

    # initialize the sequence accumulator
    sequence = ''

    # process all the rest of the lines in the file
    for line in infile:
        line = line.rstrip()

        # ignore blank lines
        if line != '':

            # if it's a header line, finish the previous sequence
            # and start a new one
            if line[0] == '>':
                resultList.append([label, headerLine, sequence])

                label = parseHeaderLine(line)
                headerLine = line.rstrip()
                sequence = ''
            
            # if we're here, we must be in letters of the sequence
            else:
                sequence += line
            
    # we're done, so clean up, terminate the last sequence, and return
    infile.close()
    resultList.append([label, headerLine, sequence])
    return resultList

main()




"#DNAalignment" 
