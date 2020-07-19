import pysam

def do(inbam, outbam, thresh):
    bamin = pysam.AlignmentFile(inbam, 'rb')
    bamout = pysam.AlignmentFile(outbam, 'wb', template=bamin)
    for read in bamin:
        if read.get_tag('ZW') > thresh:
            bamout.write(read)
    bamin.close()
    bamout.close()
    