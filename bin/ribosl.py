import re, sys, random, os, math, pickle, time, scipy
from scipy import stats
from math import isnan
from multiprocessing import Pool
from functools import partial

Bin = sys.path[0]


def ReadFasta(inFa, inDict):
    rowSize = 0
    seqAll = []
    with open(inFa, 'r') as f:
        seqAll = f.read().split('>')
    for line in seqAll:
        if len(line) == 0: continue
        Line = line.split('\n')
        idLine = Line.pop(0)
        seqLine = ''.join(Line)
        cchr = idLine.replace('>', '').split()[0]
        inDict[cchr] = seqLine


def RC(seq, type='DNA'):
    """Return reverse and complementary DNA or RNA sequences"""
    seq = seq[::-1].upper()
    for pair in ['At', 'Ta', 'Cg', 'Gc']:
        seq = seq.replace(pair[0], pair[1])
    return seq.upper() if type == 'DNA' else seq.upper().replace('T', 'U')


def Entropy3(lst):
    """Calculating entropy of a list using log3"""
    entropy = 0
    for p in lst:
        p = max(p, 0.0001)
        entropy -= p * math.log(p) / math.log(3)
    return entropy


def IndexE(seq, codon):
    """return all the positions of codon (eg. AUG, NUG)"""
    OutArr = []
    posHit = 0
    posStart = -1
    while posHit >= 0:
        posHit = seq.find(codon, posStart + 1)
        posStart = posHit
        OutArr.append(posHit)
    return OutArr


def Zscore(lst):
    """return Z-score of each number in the list"""
    Zscore = []
    for number in lst:
        z = 0 if scipy.std(lst) == 0 else (number - scipy.mean(lst)) / scipy.std(lst)
        Zscore.append(z)
    return Zscore


def StartExtract(inGff, pre):
    """Extract positions of start and stop codon from GTF"""
    fileIn = open(inGff, 'r')
    GFF = {}
    for eachLine in fileIn:
        Line = eachLine.rstrip().split()
        if Line[0].startswith('#'): continue
        if 'CDS' not in Line[2]: continue
        transId = re.findall(r'transcript_id \"([^\"]+)\"', eachLine)[0]
        cchr, start, end, strand = Line[0], int(Line[3]), int(Line[4]), Line[6]
        start -= 1  # convert to 0-base
        if transId not in GFF: GFF[transId] = []
        GFF[transId].append([cchr, start, end, strand])
    fileIn.close()

    startFile = open(pre + '.start.tmp', 'w')
    stopFile = open(pre + '.stop.tmp', 'w')
    for transId in GFF:
        startFrom, startTo, stopFrom, stopTo = 0, 0, 0, 0
        cchr, strand = GFF[transId][0][0], GFF[transId][0][3]
        ExonSort = sorted(GFF[transId], key=lambda exon: exon[1]) if strand == '+' else sorted(GFF[transId],
                                                                                               key=lambda exon: exon[1],
                                                                                               reverse=True)
        if strand == '+':
            startFrom = ExonSort[0][1]
            startTo = startFrom + 3
            stopTo = ExonSort[-1][2]
            stopFrom = stopTo - 3
        else:
            startTo = ExonSort[0][2]
            startFrom = startTo - 3
            stopFrom = ExonSort[-1][1]
            stopTo = stopFrom + 3
        print(cchr, str(startFrom), str(startTo), transId, '.', strand, sep='\t', file=startFile)
        print(cchr, str(stopFrom), str(stopTo), transId, '.', strand, sep='\t', file=stopFile)
    startFile.close()
    stopFile.close()
    os.system('sort -k1,1 -k2,2n ' + pre + '.start.tmp > ' + pre + '.start.bed')
    os.system('sort -k1,1 -k2,2n ' + pre + '.stop.tmp > ' + pre + '.stop.bed')
    os.system('rm *.tmp')


def OffSetExtract(bamFile, startCodon, pre, logout, pathBedtools='', pathSamtools=''):
    """ Calculate offset for RPFs in each size and return the entropy of frames for size with most abundance"""
    if not pathBedtools:
        out = os.system('bamToBed -bed12 -i ' + bamFile + ' > ' + pre + '.bed 2>>' + logout)
        if out != 0:
            sys.stdout.write(
                "ERROR: bedtools is not an excutable command, please provide the path to bedtools using the option of --bedtools")
            sys.exit()
    else:
        os.system(pathBedtools + '/bamToBed -bed12 -i ' + bamFile + ' > ' + pre + '.bed 2>>' + logout)

    if not pathSamtools:
        out = os.system(
            'samtools view -h -s 1.03 ' + bamFile + ' 2>>' + logout + '|samtools view -bS - > ' + pre + '.sub.bam 2>>' + logout)
        if out != 0:
            sys.stdout.write(
                "ERROR: samtools is not an excutable command, please provide the path to bedtools using the option of --samtools")
            sys.exit()
    else:
        os.system(
            pathSamtools + '/samtools view -h -s 1.03 ' + bamFile + ' 2>>' + logout + '|' + pathSamtools + '/samtools view -bS - > ' + pre + '.sub.bam 2>>' + logout)

    if not pathBedtools:
        os.system('bamToBed -bed12 -i ' + pre + '.sub.bam > ' + pre + '.sub.bed 2>>' + logout)
        os.system(
            'windowBed -w 100 -sm -b ' + pre + '.sub.bed -a ' + startCodon + ' 2>>' + logout + '|cut -f 7-18|sort -k1,1 -k2,2g|closestBed -s -a stdin -b ' + startCodon + ' > ' + pre + ' 2>>' + logout)
    else:
        os.system(pathBedtools + '/bamToBed -bed12 -i ' + pre + '.sub.bam > ' + pre + '.sub.bed 2>>' + logout)
        os.system(
            pathBedtools + '/windowBed -w 100 -sm -b ' + pre + '.sub.bed -a ' + startCodon + ' 2>>' + logout + '|cut -f 7-18|sort -k1,1 -k2,2g|' + pathBedtools + '/closestBed -s -a stdin -b ' + startCodon + ' > ' + pre + ' 2>>' + logout)
    os.system('rm *sub.bam *sub.bed')

    Off, Len = {}, {}
    fileIn = open(pre, 'r')
    for line in fileIn:
        Line = line.rstrip().split()
        if ',' in Line[10]: continue
        strand, size = Line[5], int(Line[10])
        dist = int(Line[1]) - int(Line[13]) if strand == '+' else int(Line[14]) - int(Line[2])
        if size not in Off: Off[size] = {}
        if dist not in Off[size]: Off[size][dist] = 0
        if size not in Len: Len[size] = 0
        Off[size][dist] += 1
        Len[size] += 1
    fileIn.close()
    R = open(pre + '.R.txt', 'w')
    L = open(pre + '.L.txt', 'w')
    print('ReadLen\t' + "\t".join(list(map(str, list(range(0, 51))))), file=R)
    print('ReadLen\t' + "\t".join(list(map(str, list(range(-20, 0))))), file=L)
    for size in sorted(Len):
        outR, outL = str(size), str(size)
        for dist in range(-20, 51):
            count = Off[size][dist] if dist in Off[size] else 0
            count = Off[size][dist] if dist in Off[size] else 0
            if dist >= 0:
                outR += '\t' + str(count)
            else:
                outL += '\t' + str(count)
        print(outR, file=R)
        print(outL, file=L)
    R.close()
    L.close()

    Entropy = {}
    fileIn = open(pre + '.L.txt', 'r')
    Title = list(map(int, fileIn.readline().rstrip().split('\t')[1:]))
    fileout = open(pre + '_OffSet.txt', 'w')
    for line in fileIn:
        Line = list(map(int, line.rstrip().split()))
        size = Line.pop(0)
        offset, offsetP = [], []
        maxDep = max(Line)
        for i in range(len(Line) - 2):
            count = Line[i]
            if count <= 0.8 * maxDep or size + Title[i] <= 0: continue
            if sum(Line[i:i + 3]) == 0: continue
            for pos in [0, 1, 2]:
                offset.append(abs(Title[i + pos]))
                offsetP.append(Line[i + pos] / sum(Line[i:i + 3]))
                print(size, offset[pos], offsetP[pos], sep='\t', file=fileout)
            Entropy[size] = Entropy3(offsetP)
            break
    fileIn.close()
    fileout.close()

    with open('entropy.pickle', 'wb') as f, open('SizeDis.pickle', 'wb') as ff:
        pickle.dump(Entropy, f)
        pickle.dump(Len, ff)


def read_bam(bam_file, sl, soft_num, filetype):
    soft_num = int(soft_num)
    if filetype == 'bam':
        os.system(f"samtools view -h {bam_file} | less -S > test.txt")
    else:
        os.system(f"cp {bam_file} ./test.txt")
    SLTrans = {}
    with open('test.txt', 'r') as f:
        whole_db = f.readlines()
        for single in whole_db:
            if single.startswith(('@', '#')): continue
            single = single.split()
            seq_mod = single[5]
            CigarValues = re.split("[A-Z]", seq_mod)
            slPosition, softclip = 0, 0
            slTail = ''
            if seq_mod.endswith('S'):
                softclip = int(CigarValues[-2])
                if softclip < soft_num: continue
                slPosition = int(single[3]) + sum(map(int, list(CigarValues[:-2])))
                slTail = RC(single[9][-1 * softclip:])
            elif 'S' in seq_mod:
                softclip = int(CigarValues[0])
                if softclip < soft_num: continue
                slPosition = int(single[3])
                slTail = single[9][:softclip]
            else:
                continue

            for slSeq in sl.split(','):  # multiple SL sequences
                if slSeq.endswith(slTail):
                    coordFull = (single[2], slPosition)
                elif slTail in slSeq:
                    tailExtend = len(slSeq.split(slTail)[-1])
                    # sys.stdout.write(slSeq + '\t\t' + slTail + str(tailExtend) +'\n')
                    coordFull = (single[2], slPosition - tailExtend) if seq_mod.endswith('S') else (
                    single[2], slPosition + tailExtend)
                else:
                    continue
                if coordFull not in SLTrans: SLTrans[coordFull] = {}
                if slSeq not in SLTrans[coordFull]: SLTrans[coordFull][slSeq] = 0
                SLTrans[coordFull][slSeq] += 1

    os.system('rm -f test.txt')
    return (SLTrans)


def codonusage(transId, GFF, SEQ):
    Codon = {}
    startFrom, startTo, stopFrom, stopTo = 0, 0, 0, 0
    cchr, strand = GFF[transId][0][0], GFF[transId][0][3]
    ExonSort = sorted(GFF[transId], key=lambda exon: exon[1]) if strand == '+' else sorted(GFF[transId],
                                                                                           key=lambda exon: exon[1],
                                                                                           reverse=True)
    transIdNew = '\t'.join([transId, cchr, strand])
    cds = ''
    for cchr, start, end, strand in ExonSort:
        exon = SEQ[cchr][start - 1:end] if strand == '+' else RC(SEQ[cchr][start - 1:end])
        exon = exon.replace('T', 'U')
        cds += exon
    for i in range(0, len(cds), 3):
        codon = cds[i:i + 3]
        if len(re.findall(re.compile('A|T|U|C|G'), codon)) != 3: continue
        if codon not in Codon: Codon[codon] = 0
        Codon[codon] += 1
    return (Codon)


def UsageExtract(genome, inGtf, nCores, pre):
    """Calculate the genome-wide codon usage for each codon"""
    SEQ, GFF, Codon, Usage = {}, {}, {}, {}
    codonSum = 0
    ReadFasta(genome, SEQ)
    with open(pre + '_ref.pickle', 'wb') as f:
        pickle.dump(SEQ, f)

    gffFile = open(inGtf, 'r')
    for eachLine in gffFile:
        Line = eachLine.rstrip().split()
        if Line[0].startswith('#'): continue
        if 'CDS' not in Line[2]: continue
        transId = re.findall(r'transcript_id \"([^\"]+)\"', eachLine)[0]
        cchr, start, end, strand = Line[0], int(Line[3]), int(Line[4]), Line[6]
        if transId not in GFF: GFF[transId] = []
        GFF[transId].append([cchr, start, end, strand])
    gffFile.close()

    for transId in GFF:
        startFrom, startTo, stopFrom, stopTo = 0, 0, 0, 0
        cchr, strand = GFF[transId][0][0], GFF[transId][0][3]
        ExonSort = sorted(GFF[transId], key=lambda exon: exon[1]) if strand == '+' else sorted(GFF[transId],
                                                                                               key=lambda exon: exon[1],
                                                                                               reverse=True)
        transIdNew = '\t'.join([transId, cchr, strand])
        cds = ''
        for cchr, start, end, strand in ExonSort:
            exon = SEQ[cchr][start - 1:end] if strand == '+' else RC(SEQ[cchr][start - 1:end])
            exon = exon.replace('T', 'U')
            cds += exon
        for i in range(0, len(cds), 3):
            codon = cds[i:i + 3]
            if len(re.findall(re.compile('A|T|U|C|G'), codon)) != 3: continue
            if codon not in Codon: Codon[codon] = 0
            Codon[codon] += 1
            codonSum += 1
    for codon in Codon:
        usage = Codon[codon] / codonSum
        codon.replace('T', 'U')
        Usage[codon] = usage
    with open(pre + '_codonUsage.pickle', 'wb') as f:
        pickle.dump(Usage, f)


def PsiteAllocate(inBed, inOffset, inFrq, pre):
    """Allocate the RPFs to Psites they are representing according to their offsets and coresponding probabilities."""
    Entropy, SizeDis = {}, {}
    SizeFrqH, SizeFrqL, SizeSelect = [], [], []
    numUsedReadH, numUsedReadL, numUsedRead, numDiscardRead, numTotalRead = 0, 0, 0, 0, 0
    with open('entropy.pickle', 'rb') as f, open('SizeDis.pickle', 'rb') as ff:
        Entropy = pickle.load(f)
        SizeDis = pickle.load(ff)

    fileIn = open(inFrq, 'r')
    for line in fileIn:
        if line.startswith('#'): continue
        size, logP, frq = line.rstrip().split()
        size, logP, frq = int(size), float(logP), float(frq)
        if (frq < 0.33 or frq > 0.34 or logP < 1.3) and (size in Entropy): del Entropy[size]
    fileIn.close()

    sumAmount, sumEntropy, meanEntropy = 0, 0, 0
    for size, entropy in sorted(Entropy.items(), key=lambda val: val[1])[0:4]:  # Top 5 periodic RPFs
        amount = SizeDis[size]
        sumAmount += amount
        sumEntropy += entropy * amount
        SizeSelect.append(size)
    if sumAmount == 0:
        sys.stdout.write("\tThere is no periodical RPFs can be used for prediction. QUIT!" + '\n')
        quit()
    meanEntropy = sumEntropy / sumAmount

    Prob = {}
    fileIn = open(inOffset, 'r')
    for line in fileIn:
        size, offset, prob = line.rstrip().split()  # prob is short for probability
        size = int(size)
        offset = int(offset)
        prob = float(prob)
        offsetAnti = size - offset - 1
        if size not in Prob: Prob[size] = {}
        Prob[size][('+', offset)] = prob
        Prob[size][('-', offsetAnti)] = prob
    fileIn.close()

    fileIn = open(inBed, 'r')
    for line in fileIn:
        Line = line.rstrip().split()
        size = sum(list(map(int, Line[10].split(','))))
        if size in SizeSelect: numUsedRead += 1
        numTotalRead += 1
    fileIn.close()

    numDiscardRead = numTotalRead - numUsedRead

    PprobSum = {}
    fileIn = open(inBed, 'r')
    for line in fileIn:
        Line = line.rstrip().split()
        Sizes = list(map(int, Line[10].rstrip().split(',')))
        Starts = list(map(int, Line[11].rstrip().split(',')))
        size = sum(Sizes)
        cchr, strand = Line[0], Line[5]
        if size not in SizeSelect: continue
        if size not in Prob: continue
        for strandd, offset in sorted(Prob[size]):
            if strandd != strand: continue
            dist, psite, pprob = 0, 0, 0
            for i in range(len(Sizes)):
                if offset <= Sizes[i] + dist:
                    psite = int(Line[1]) + Starts[i] + offset - dist
                    pprob = pprob if pprob > Prob[size][(strand, offset)] else Prob[size][(strand, offset)]
                else:
                    dist += Sizes[i]
            psite += 1
            if (cchr, strand, psite) not in PprobSum: PprobSum[(cchr, strand, psite)] = 0
            PprobSum[(cchr, strand, psite)] += pprob
    fileIn.close()

    with open(pre + '_dep.pickle', 'wb') as f:
        pickle.dump(PprobSum, f)

    fileOut = open('DataSumForPrediction.txt', 'w')
    for size in sorted(SizeSelect):
        print('RPFs of ' + str(size) + ' nt were used to predict ORFs', file=fileOut)

    wz = 1 - float(meanEntropy)
    wu = float(meanEntropy)
    print('Total amount of used RPFs: ' + str(numUsedRead), file=fileOut)
    print('Total amount of discarded RPFs: ' + str(numDiscardRead), file=fileOut)
    print('Percentage of used RPFs: ' + "%.4f" % (100 * numUsedRead / numTotalRead), file=fileOut)
    print('The entropy is: ' + str(meanEntropy), file=fileOut)
    print('The weight of RPFs is: ' + str(wz), file=fileOut)
    print('The weight of codon usage is: ' + str(wu), file=fileOut)
    fileOut.close()
    return (meanEntropy)


def finder(PreOrf, pcov=0.00001, startCodon='AUG', entropy=1):
    """ predict ORFs on transcripts"""
    idLine, seqLine, coordLine, depLine, usageLine = PreOrf.split('\n')
    wz = 0.5 * (1 - float(entropy))  # weight of RPFs
    wu = 0.5 * float(entropy)  # weight of codon usage
    orfOut = ''
    Coord = coordLine.split()
    Dep = list(map(float, depLine.split()))
    Usage = list(map(float, usageLine.split()))
    ID, cchr, strand, slSeq = idLine.split('\t')
    startCodon = startCodon.replace("T", "U")
    StartCodon = startCodon.split(",")
    StopCodon = ["UGA", "UAG", "UAA"]
    PosNTG, PosTGA = [], []
    PosATG = IndexE(seqLine, 'AUG')

    for codon in StartCodon:
        if codon == 'AUG': continue
        PosHit = IndexE(seqLine, codon)
        PosNTG.extend(PosHit)

    for codon in StopCodon:
        PosHit = IndexE(seqLine, codon)
        PosTGA.extend(PosHit)

    PosATG.sort()
    PosNTG.sort()
    PosATG.extend(PosNTG)  # preORF initiated with AUG was calculated with higher priority
    PosTGA.sort()

    PosFrameATG, PosFrameTGA = [[], [], []], [[], [], []]
    for atg in PosATG: PosFrameATG[atg % 3].append(atg)
    for tga in PosTGA: PosFrameTGA[tga % 3].append(tga)

    for frameTGA in (0, 1, 2):
        for frameATG in (0, 1, 2):
            if frameATG != frameTGA: continue
            frameIn = frameATG
            for i in range(len(PosFrameTGA[frameTGA])):
                tga = PosFrameTGA[frameATG][i]
                if tga < 0: continue
                PredictedStop = []

                for j in range(len(PosFrameATG[frameATG])):
                    atg = PosFrameATG[frameATG][j]
                    if atg < 0: continue
                    if tga < atg + 45: continue  # minimum length of ORF is 45 bp
                    if (tga - atg) % 3 != 0: continue
                    if i > 0 and atg < PosFrameTGA[frameTGA][i - 1] and (
                            PosFrameTGA[frameTGA][i - 1] - atg) % 3 == 0: continue
                    if '-' not in slSeq and atg > len(slSeq) - 1: continue

                    orfSeq = seqLine[atg:tga + 3]
                    orfCoord = Coord[atg:tga + 3]
                    coordFrom = Coord[atg] if Coord[atg] < Coord[tga + 2] else Coord[tga + 2]
                    coordTo = Coord[atg] if Coord[atg] > Coord[tga + 2] else Coord[tga + 2]
                    orfId = ":".join([ID, cchr, strand, str(coordFrom), str(coordTo)])
                    '''
                    if '_' in ID and atg > len(SL):
                        pass
                    else:
                    '''
                    DepOrfAllsites = Dep[atg:tga + 3]
                    if sum(DepOrfAllsites) == 0: continue
                    DepZ = Zscore(DepOrfAllsites)
                    ZOrf = [[], [], []]  # Z-score at each position belonging to frame 0, 1 and 2
                    UOrf = [[], [], []]

                    sumDep, pos = 0, 0
                    for k in range(atg, tga):
                        pos += 1
                        dep = Dep[k]
                        usage = Usage[k]
                        zscore = DepZ[pos - 1]
                        if dep > 0 and (pos - 1) % 3 == 0: sumDep += 1  # total depth at P-sites
                        ZOrf[(pos - 1) % 3].append(zscore)
                        UOrf[(pos - 1) % 3].append(usage)

                    if (sumDep / pos < float(
                            pcov)): continue  # quit if the coverage of Psites is smaller then pcov [default: 0]
                    tZ_1, pZ_1 = stats.ttest_ind(ZOrf[0], ZOrf[1], nan_policy='omit');
                    pZ_1 = 1 if tZ_1 < 0 or isnan(pZ_1) else 0.5 * pZ_1 + 1e-11
                    tZ_2, pZ_2 = stats.ttest_ind(ZOrf[0], ZOrf[2], nan_policy='omit');
                    pZ_2 = 1 if tZ_2 < 0 or isnan(pZ_2) else 0.5 * pZ_2 + 1e-11
                    tU_1, pU_1 = stats.ttest_ind(UOrf[0], UOrf[1], nan_policy='omit');
                    pU_1 = 1 if tU_1 < 0 or isnan(pU_1) else 0.5 * pU_1 + 1e-11
                    tU_2, pU_2 = stats.ttest_ind(UOrf[0], UOrf[2], nan_policy='omit');
                    pU_2 = 1 if tU_2 < 0 or isnan(pU_2) else 0.5 * pU_2 + 1e-11
                    freedom = getFreedom([pZ_1, pZ_2, pU_1, pU_2], [wz, wz, wu, wu])
                    pVal = 1 - stats.chi2.cdf(-2 * 0.5 * freedom * (
                            wz * math.log(pZ_1) + wz * math.log(pZ_2) + wu * math.log(pU_1) + wu * math.log(pU_2)),
                                              freedom)
                    if pVal > 0.001: continue
                    if tga in PredictedStop: continue
                    orfIdLine = '\t'.join(['>' + orfId, cchr, strand, str(pVal), ID,
                                           str(pZ_1) + ':' + str(pZ_2) + ':' + str(pU_1) + ':' + str(pU_2)])
                    orfCoordLine = ' '.join(list(map(str, orfCoord)))

                    orfdeep = ' '.join(list(map(str, Dep[atg:tga + 3])))
                    orfusage = ' '.join(list(map(str, Usage[atg:tga + 3])))

                    orfOut += orfIdLine + '\n' + orfSeq + '\n' + orfCoordLine + '\n' + orfdeep + '\n' + orfusage + '\n'
                    PredictedStop.append(tga)
                    break
    return (orfOut)


def getFreedom(Pval, Weight):
    Rho = {}
    Sm = []
    sumSm, qt = 0, 0
    for i in range(len(Pval)):
        si = -2 * math.log(Pval[i])
        Sm.append(si)
        sumSm += si
    smMean = sumSm / len(Pval)

    for i in range(len(Sm)):
        si = Sm[i]
        qt += (si - smMean) ** 2 / (len(Sm) - 1)

    for i in range(len(Pval)):
        for j in range(len(Pval)):
            if i == j: continue
            Rho[(i, j)] = -2.167 + (10.028 - 4 * qt / 3) ** 0.5 if 4 * qt / 3 < 10.028 else 0
            Rho[(i, j)] = 0 if Rho[(i, j)] < 0 else Rho[(i, j)]
            Rho[(i, j)] = 1 if Rho[(i, j)] > 1 else Rho[(i, j)]

    var = 0
    for i in range(len(Weight)):
        var += 4 * Weight[i] ** 2
        for j in range(len(Weight)):
            if i == j: continue
            rho = Rho[(i, j)]
            var += Weight[i] * Weight[j] * (3.25 * rho + 0.75 * rho ** 2)
    freedom = 8 / var
    return (freedom)


def fdr(Orfs):
    Orfs = list(filter(len, Orfs))
    OrfsOut = []
    OrfP = {}
    for orfs in Orfs:
        Line = orfs.split('\n')
        for i in range(0, len(Line) - 1, 5):
            idLine = Line[i].rstrip()
            seqLine = Line[i + 1].rstrip()
            coordLine = Line[i + 2].rstrip()

            depline = Line[i + 3].rstrip()
            usageline = Line[i + 4].rstrip()

            pval = float(idLine.split('\t')[3])
            orfInfo = idLine + '\n' + seqLine + '\n' + coordLine + '\n' + depline + '\n' + usageline + '\n'
            OrfP[orfInfo] = pval

    k = len(OrfP)
    qvalLast = 1
    for orf, pval in sorted(OrfP.items(), key=lambda val: val[1], reverse=True):
        qval = pval * len(OrfP) / k
        k -= 1
        qval = min(qval, qvalLast)
        qvalLast = qval
        if qval > 0.0001: continue  ## control fdr at 0.01%
        OrfsOut.append(orf)
    return (OrfsOut)


def TransSplicing(transcript, Coord, i, sl):
    transcript = list(transcript)
    transcript = "".join(transcript[i:])
    Coord = Coord[i:]
    sl_transcropt = sl + transcript
    sl_transcropt = sl_transcropt.replace('T', 'U')
    Coord.reverse()
    for i in range(len(list(sl))):
        Coord.append(0)
    Coord.reverse()
    return sl_transcropt, Coord


def OrfFinder(entropy, startCodon, inRef, inGtf, pcov, nCores, SLTrans, transOnly, pre):
    Ref, Exon, USAGE, Pdep, PreOrf = {}, {}, {}, {}, []
    with open(pre + '_codonUsage.pickle', 'rb') as f:
        USAGE = pickle.load(f)
    with open(pre + '_dep.pickle', 'rb') as f:
        Pdep = pickle.load(f)
    with open(pre + '_ref.pickle', 'rb') as f:
        Ref = pickle.load(f)

    fileIn = open(inGtf, 'r')
    for line in fileIn:
        if "#" not in line:
            Line = line.rstrip().split()
            if 'exon' not in Line[2]: continue
            cchr, exonFrom, exonTo, strand = Line[0], int(Line[3]), int(Line[4]), Line[6]
            transId = re.findall(r'transcript_id \"([^\"]+)\"', line)[0]
            if transId not in Exon: Exon[transId] = []
            Exon[transId].append([cchr, strand, exonFrom, exonTo])
    fileIn.close()

    for transId in Exon:
        Rpt = []
        transcript = ''
        Coord = []
        cchr, strand = Exon[transId][0][0:2]
        ExonSort = sorted(Exon[transId], key=lambda exon: exon[2]) if strand == '+' else sorted(Exon[transId],
                                                                                                key=lambda exon: exon[
                                                                                                    2], reverse=True)
        for exon in ExonSort:
            if exon in Rpt: continue
            Rpt.append(exon)
            exonFromE, exonToE = exon[2:4]
            seqExon = Ref[cchr][exonFromE - 1:exonToE]
            coordExon = list(range(exonFromE, exonToE + 1)) if strand == '+' else list(
                range(exonToE, exonFromE - 1, -1))
            seqExon = RC(seqExon, 'RNA') if strand == '-' else seqExon.replace('T', 'U')
            Coord.extend(coordExon)
            transcript += seqExon

        ### SL trans-splicing ####
        SLtranscript = [transcript]
        SLcoord = [Coord]
        SLSeq = ['-']
        for i in range(len(transcript)):
            coord = Coord[i]
            coordFull = (cchr, coord)
            if coordFull in SLTrans:
                for slSeq in SLTrans[coordFull]:
                    if SLTrans[coordFull][slSeq] < 1: continue
                    sltranscript, slcoord = TransSplicing(transcript, Coord, i, slSeq)
                    # print(f'\n {coordFull} \n {sltranscript}')
                    SLtranscript.append(sltranscript)
                    SLcoord.append(slcoord)
                    SLSeq.append(slSeq)
        if len(SLtranscript) == transOnly: continue
        for k in range(transOnly, len(SLtranscript)):  # TransOnly: 0/1，1: only transORF, 0: all ORFs
            Dep, Usage = [], []
            slSeq = SLSeq[k]
            for i in range(len(SLtranscript[k])):
                coord = SLcoord[k][i]
                coordFull = (cchr, strand, coord)
                dep = Pdep[coordFull] if coordFull in Pdep else 0
                codon = transcript[i:i + 3]
                usage = USAGE[codon] if codon in USAGE else 0
                Dep.append(dep)
                Usage.append(usage)
            if k > 0:
                transId_sl = transId + "_sl_" + str(k)
                idLine = transId_sl + '\t' + cchr + '\t' + strand + '\t' + slSeq
            else:
                idLine = transId + '\t' + cchr + '\t' + strand + '\t' + slSeq

            if sum(Dep) == 0: continue
            seqLine = SLtranscript[k]
            coordLine = ' '.join(list(map(str, SLcoord[k])))
            depLine = ' '.join(list(map(str, Dep)))
            usageLine = ' '.join(list(map(str, Usage)))
            PreOrf.append(idLine + '\n' + seqLine + '\n' + coordLine + '\n' + depLine + '\n' + usageLine)
            #print(idLine + '\n' + seqLine + '\n' + coordLine + '\n' + depLine + '\n' + usageLine)

    if (nCores > 1):
        with Pool(processes=nCores) as pool:
            Orfs = pool.map(partial(finder, pcov=pcov, startCodon=startCodon, entropy=entropy), PreOrf)
    else:
        Orfs = map(partial(finder, pcov=pcov, startCodon=startCodon, entropy=entropy), PreOrf)

    if transOnly == 0:
        Orfs = list(filter(len,Orfs))
    else:
        Orfs = fdr(Orfs)
    with open(pre + '_orf.fa', 'w') as f:
        f.write(''.join(filter(len, Orfs)))


def gffMaker(inGtf, startBed, stopBed, inCds, fileOutGff, outFa, outTable, outSum):
    CODE = {
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',  # Alanine
        'UGC': 'C', 'UGU': 'C',  # Cysteine
        'GAC': 'D', 'GAU': 'D',  # Aspartic Acid
        'GAA': 'E', 'GAG': 'E',  # Glutamic Acid
        'UUC': 'F', 'UUU': 'F',  # Phenylalanine
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G',  # Glycine
        'CAC': 'H', 'CAU': 'H',  # Histidine
        'AUA': 'I', 'AUC': 'I', 'AUU': 'I',  # Isoleucine
        'AAA': 'K', 'AAG': 'K',  # Lysine
        'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L', 'UUA': 'L', 'UUG': 'L',  # Leucine
        'AUG': 'M',  # Methionine
        'AAC': 'N', 'AAU': 'N',  # Asparagine
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',  # Proline
        'CAA': 'Q', 'CAG': 'Q',  # Glutamine
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R', 'AGA': 'R', 'AGG': 'R',  # Arginine
        'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S', 'AGC': 'S', 'AGU': 'S',  # Serine
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',  # Threonine
        'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V',  # Valine
        'UGG': 'W',  # Tryptophan
        'UAC': 'Y', 'UAU': 'Y',  # Tyrosine
        'UAA': 'U', 'UAG': 'U', 'UGA': 'U'  # Stop
    }

    EleAll, Ele, OrfClass, Anno = {}, {}, {}, {}
    ## Load the position of annotated positions of start codon
    fileIn = open(startBed, 'r')
    for line in fileIn:
        cchr, start, end, transId, score, strand = line.rstrip().split()
        start, end = int(start), int(end)
        start += 1  # convert to 1-based coordinates
        if transId not in Anno: Anno[transId] = {}
        Anno[transId]['start'] = [cchr, strand, start] if strand == '+' else [cchr, strand, end]
    fileIn.close()

    ## Load the position of annotated positions of stop codon
    fileIn = open(stopBed, 'r')
    for line in fileIn:
        cchr, start, end, transId, score, strand = line.rstrip().split()
        start, end = int(start), int(end)
        start += 1  ## convert to 1-based coordinates
        if transId not in Anno: Anno[transId] = {}
        Anno[transId]['stop'] = [cchr, strand, end] if strand == '+' else [cchr, strand, start]
    fileIn.close()

    ## Load genomic elements for each transcript
    fileIn = open(inGtf, 'r')
    for line in fileIn:
        try:
            if 'transcript_id' not in line: continue
            Line = line.rstrip().split('\t')
            if re.findall(re.compile('UTR|CDS|exon|RNA|gen|Gen'), Line[2]):
                transId = re.findall(r'transcript_id \"([^\"]+)\"', line)[0]
                if (transId, Line[2]) not in Ele: Ele[(transId, Line[2])] = []
                Ele[(transId, Line[2])].extend([Line[3], Line[4]])
        except:
            pass
    fileIn.close()

    for transId, ele in Ele:
        geneId = transId.split('.')[0]
        if transId not in EleAll: EleAll[transId] = ''
        EleAll[transId] += ele
        if geneId not in EleAll: EleAll[geneId] = ''
        EleAll[geneId] += ele

    fileIn = open(inCds, 'r')
    fileOut = open(fileOutGff, 'w')
    line = '1'
    while line:
        line = fileIn.readline().rstrip()
        if not line.startswith('>'): continue
        idLine = line
        seqLine = fileIn.readline().rstrip()
        coordLine = fileIn.readline().rstrip()
        Coord = list(map(int, coordLine.split()))
        orfId, cchr, strand, orfPval, transId, info = idLine.replace('>', '').split('\t')
        Start, End = [], []
        ii = 0
        for i in range(len(Coord)):
            if Coord[i] == 0:
                ii += 1
            else:
                if i - ii == 0:
                    Start.append(Coord[i])
                elif i == len(Coord) - 1:
                    End.append(Coord[i])
                elif strand == '+' and Coord[i] + 1 != Coord[i + 1]:
                    Start.append(Coord[i + 1])
                    End.append(Coord[i])
                elif strand == '-' and Coord[i] - 1 != Coord[i + 1]:
                    Start.append(Coord[i + 1])
                    End.append(Coord[i])
        '''
        修改0（只修改了第一列，没有修改剩余列）
        # 修改gff，防止起始密码子坐标为0
        for start_num in range(len(Start)):
            if int(Start[start_num]) != 0:
                startOrf = Start[start_num]
        for end_num in range(1, len(Start)):
            if int(End[-1*end_num]) != 0:
                endOrf = End[-1*end_num]
        
        # 修改gff，防止起始密码子坐标为0（将0全部去除）
        sa = Start.count(0)
        en = End.count(0)
        for start_num in range(0, sa):
            Start.remove(0)
        for end_num in range(0, en):
            End.remove(0)
        '''

        startOrf, endOrf = int(Start[0]), int(End[-1])

        # print(startOrf)
        # print(Start)
        orfClass = ''
        geneId = transId.split('.')[0]

        if transId in Anno:
            annoStart = Anno[transId]['start'][-1]
            annoStop = Anno[transId]['stop'][-1]
            if strand == '+':
                if endOrf == annoStop or endOrf == annoStop + 3:  ## Stop codon is not included in CDSs in some GTFs, such as human
                    if startOrf == annoStart: orfClass = 'annotated_ORF'
                    if startOrf > annoStart: orfClass = 'truncated_ORF'
                    if startOrf < annoStart: orfClass = 'extended_ORF'
                elif startOrf > annoStop:
                    orfClass = 'dORF'  ## downstream ORF
                elif endOrf > annoStop:
                    orfClass = 'odORF' if startOrf > annoStart else 'novel_ORF'  ## odORF: overlapped dORF
                elif endOrf <= annoStart:
                    orfClass = 'uORF'
                elif endOrf < annoStop:
                    orfClass = 'ouORF' if startOrf < annoStart else 'internal_ORF'  ## ouORF: overlapped uORF
            else:
                if endOrf == annoStop or endOrf == annoStop - 3:
                    if startOrf == annoStart: orfClass = 'annotated_ORF'
                    if startOrf < annoStart: orfClass = 'truncated_ORF'
                    if startOrf > annoStart: orfClass = 'extended_ORF'
                elif startOrf < annoStop:
                    orfClass = 'dORF'
                elif endOrf < annoStop:
                    orfClass = 'odORF' if startOrf < annoStart else 'novel_ORF'
                elif endOrf >= annoStart:
                    orfClass = 'uORF'
                elif endOrf > annoStop:
                    orfClass = 'ouORF' if startOrf > annoStart else 'internal_ORF'
        elif geneId in EleAll:
            if 'CDS' in EleAll[geneId]:
                orfClass = 'novel_ORF'
            elif 'TE' in EleAll[geneId]:
                orfClass = 'teORF'
            elif 'pseudo' in EleAll[geneId]:
                orfClass = 'pORF'
            else:
                orfClass = 'ncsORF'
        else:
            orfClass = 'intergenic_ORF'
        if '_sl_' in transId:
            orfClass = 'sl_transORF'
        OrfClass[orfId] = orfClass
        print(cchr + '\tslORFfinder\t' + orfClass + '\t' + '\t'.join(
            list(map(str, sorted([startOrf, endOrf])))) + '\t' + orfPval + '\t' + strand + '\t.\tID=' + orfId,
              file=fileOut)
        for i in range(len(Start)):
            start = Start[i] if strand == '+' else End[i]
            end = End[i] if strand == '+' else Start[i]
            print(cchr + '\tslORFfinder\texon\t' + str(start) + '\t' + str(end) + '\t.\t' + strand + '\t.\tParent=' + orfId,
                  file=fileOut)
    fileIn.close()
    fileOut.close()

    fileIn = open(inCds, 'r')
    fileOutA = open(outFa, 'w')
    fileOutT = open(outTable, 'w')
    NumClass = {}
    print('OrfId\tClass\tChr\tStrand\tPval\tCDS\tAA', file=fileOutT)
    line = '1'
    while line:
        line = fileIn.readline().rstrip()
        if not line.startswith('>'): continue
        idLine = line
        seqLine = fileIn.readline().rstrip()
        coordLine = fileIn.readline().rstrip()
        Coord = coordLine.split()
        orfId, cchr, strand, orfPval, transId, info = idLine.replace('>', '').split()
        orfClass = OrfClass[orfId]
        if orfClass not in NumClass:
            NumClass[orfClass] = 0
        else:
            NumClass[orfClass] += 1

        seqAa = ''
        for i in range(0, len(seqLine), 3):
            codon = seqLine[i:i + 3]
            aa = CODE[codon] if codon in CODE else 'X'
            seqAa += aa
        print('>' + orfId + '\n' + seqAa, file=fileOutA)
        print(orfId, orfClass, cchr, strand, orfPval, seqLine, seqAa, sep='\t', file=fileOutT)
    fileIn.close()
    fileOutA.close()
    fileOutT.close()

    fileOutS = open(outSum, 'w')
    for orfClass in NumClass:
        print(orfClass, str(NumClass[orfClass]), sep="\t", file=fileOutS)
    fileOutS.close()
