import sys
before = sys.argv[1]
after = sys.argv[2]
bef = open(before)
aft = open(after)
count = 0
for bline in bef:
    bline = bline.rstrip()
    count += 1
    aline = aft.readline().rstrip()
    print()
    blarr = bline.split('\t')
    id_ori = blarr[1]
    mid = id_ori.split(' ')[0]
    ori = id_ori.split(' ')[2].replace(')','')
    alarr = aline.split('\t')
    if mid != alarr[0]:
        raise Exception('lines don\'t match: ',mid,alarr[0])
    print('%s\t%s\t%s\tbrorig=%s\t%s\t%s\t%s' % (mid,alarr[1],alarr[2],ori,alarr[3],alarr[4],alarr[5]))
    if count == 1:
        break
bef.close()
aft.close()

