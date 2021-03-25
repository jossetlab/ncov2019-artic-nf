from __future__ import print_function
import os
import sys
import getopt

##Count the number of minor variants in a target vcf reported as major variant in a reference vcf
#v0.01

def main(argv):
    global ref
    global var
    global con
    global oup
    global mode
    global bed
    global min_depth
    global min_freq
    ref = ''
    var = ''
    con = ''
    oup = ''
    mode = ''
    bed = ''
    min_depth = 20
    min_freq = 0.01
    try:
        opts, args = getopt.getopt(argv, 'hr:c:v:o:m:b:d:f:', ['help', 'ref', 'con', 'var', 'output', 'mode', 'bed', 'min_depth', 'min_freq'])
        for opt, arg in opts:
            if opt in ('-h', '--help'):
                usage()
                sys.exit()
            elif opt in ('-r', '--ref'):
                ref = arg
            elif opt in ('-c', '--con'):
                con = arg
            elif opt in ('-v', '--var'):
                var = arg
            elif opt in ('-o', '--output'):
                oup = arg
            elif opt in ('-m', '--mode'):
                mode = arg
            elif opt in ('-b', '--bed'):
                bed = arg
            elif opt in ('-d', '--min_depth'):
                min_depth = int(arg)
            elif opt in ('-f', '--min_freq'):
                min_freq = float(arg)
        if ref == '' or con == '' or var == '':
            usage()
            sys.exit()
        if oup == '':
            oup = var.split("/")[-1].rstrip('.' + var.split('.')[-1]) + '_' + con.split("/")[-1].rstrip('.' + con.split('.')[-1]) + '.tsv'
    except getopt.GetoptError:
        usage()
        sys.exit(2)

def usage():
    print('usage: ' + sys.argv[0] + ' -h --help -r --ref [fasta] --con [vcf] --var [vcf] -o --output [tsv] -m --mode [raw,cov] -b --bed [bed] -d --min_depth [int] -f --min_freq [float]')

if __name__ == '__main__':
    main(sys.argv[1:])

def count_commented(file):
    lines = open(file, 'r').read().rstrip('\n').split('\n')
    count = 0
    for line in lines:
        if line[0] == "#":
            count += 1
    return count

flatten = lambda t: [item for sublist in t for item in sublist]

seq = [[x.replace('\r\n','\n').split('\n')[0], ''.join(x.replace('\r\n','\n').split('\n')[1:]).replace(' ','')] for x in open(ref, 'r').read().rstrip('\n').split('>')[1:]]
#cons = open(con, 'r').read().rstrip('\n').split('\n')[count_commented(con):]
#vas = open(var, 'r').read().rstrip('\n').split('\n')[count_commented(var):]
cons = open(con, 'r').read().rstrip('\n').split('\n')[1:]
vas = open(var, 'r').read().rstrip('\n').split('\n')[1:]
bga = [x.split('\t') for x in open(bed, 'r').read().replace('\r\n','\n').rstrip('\n').split('\n')]

depth = []
for i in range(len(bga)):
    depth.append([int(bga[i][3]) for x in range(int(bga[i][1]),int(bga[i][2]))])
depth = flatten(depth)

#indexs = open(con, 'r').read().split('\n')[count_commented(con)-1:count_commented(con)].split('\t')

vas_chrom, vas_pos, vas_ref, vas_alt, vas_af, cons_chrom, cons_pos, cons_ref, cons_alt, cons_af = ([] for i in range(10))
temp = []
exp = []
common = 0
expected = 0

for i in range(len(vas)):
    vas_chrom.append(vas[i].split('\t')[0])
    vas_pos.append(int(vas[i].split('\t')[1])-1)
#    vas_ref.append(vas[i].split('\t')[3])
#    vas_alt.append(vas[i].split('\t')[4])
    vas_ref.append(vas[i].split('\t')[2])
    vas_alt.append(vas[i].split('\t')[3])
#    vas_af.append(vas[i].split('\t')[7].split(';')[3].split('=')[1])
    vas_af.append(vas[i].split('\t')[10])

for i in range(len(cons)):
    cons_chrom.append(cons[i].split('\t')[0])
    cons_pos.append(int(cons[i].split('\t')[1])-1)
#    cons_ref.append(cons[i].split('\t')[3])
#    cons_alt.append(cons[i].split('\t')[4])
    cons_ref.append(cons[i].split('\t')[2])
    cons_alt.append(cons[i].split('\t')[3])
#    cons_af.append(cons[i].split('\t')[7].split(';')[3].split('=')[1])
    cons_af.append(cons[i].split('\t')[10])

for i in range(len(cons_chrom)):
    if cons_alt[i][0] == '-':
        cons_temp = cons_ref[i]
        cons_ref[i] = cons_ref[i] + cons_alt[i][1:]
        cons_alt[i] = cons_temp
    if cons_alt[i][0] == '+':
        cons_alt[i] = cons_ref[i] + cons_alt[i][1:]

for i in range(len(vas_chrom)):
    if vas_alt[i][0] == '-':
        vas_temp = vas_ref[i]
        vas_ref[i] = vas_ref[i] + vas_alt[i][1:]
        vas_alt[i] = vas_temp
    if vas_alt[i][0] == '+':
        vas_alt[i] = vas_ref[i] + vas_alt[i][1:]

for i in range(len(cons_chrom)):
    if float(cons_af[i]) >= 0.50 and float(cons_af[i]) <= 0.95:
        cons_chrom.append(cons_chrom[i])
        cons_pos.append(cons_pos[i])
        cons_ref.append(cons_alt[i])
        cons_alt.append(cons_ref[i])
        cons_af.append(1-float(cons_af[i]))

for i in range(len(vas_chrom)):
    if float(vas_af[i]) >= 0.50 and float(vas_af[i]) <= 0.95:
        vas_chrom.append(vas_chrom[i])
        vas_pos.append(vas_pos[i])
        vas_ref.append(vas_alt[i])
        vas_alt.append(vas_ref[i])
        vas_af.append(1-float(vas_af[i]))

for i in range(len(cons_chrom)):
    if depth[cons_pos[i]] >= min_depth and float(cons_af[i]) >= min_freq:
        if float(cons_af[i]) >= 0.5:
            if cons_pos[i] in vas_pos:
                vas_index = vas_pos.index(cons_pos[i])
                if cons_alt[i] == vas_alt[vas_index] and float(vas_af[vas_index]) >= 0.5:
                    print([cons_pos[i], cons_alt[i], vas_alt[vas_index], "old"])
                else:
                    expected += 1
                    exp.append([cons_pos[i], cons_alt[i], vas_alt[vas_index]])
            else:
                expected += 1
                exp.append([cons_pos[i], cons_alt[i], "ref1"])

for i in range(len(vas_chrom)):
    if depth[vas_pos[i]] >= min_depth and float(vas_af[i]) >= min_freq:
        if float(vas_af[i]) < 0.5:
            if vas_pos[i] in cons_pos:
                cons_index = cons_pos.index(vas_pos[i])
                if vas_alt[i] == cons_alt[cons_index] and float(cons_af[cons_index]) >= 0.5:
                    common += 1
                    temp.append([vas_pos[i], vas_alt[i], cons_alt[cons_index]])
            elif vas_alt[i] == seq[[x[0] for x in seq].index(vas_chrom[i])][1][vas_pos[i]:vas_pos[i]+len(vas_alt[i])] and vas_ref[i] != seq[[x[0] for x in seq].index(vas_chrom[i])][1][vas_pos[i]:vas_pos[i]+len(vas_ref[i])]:
                common += 1
                temp.append([vas_pos[i], vas_alt[i], "ref2"])
        else:
            if vas_pos[i] in cons_pos:
                cons_index = cons_pos.index(vas_pos[i])
                if vas_alt[i] != cons_alt[cons_index] and float(cons_af[cons_index]) >= 0.5:
                    expected += 1
                    exp.append([vas_pos[i], vas_alt[i], cons_alt[cons_index]])
            else:
                expected += 1
                exp.append([vas_pos[i], vas_alt[i], "ref3"])


w = open(oup, 'a+')
exp = sorted(exp, key=lambda i: i[0])
print(exp)
print (temp)
if mode == "cov":
    if expected > 0:
        w.write(con.split("/")[-1].rstrip('.' + con.split('.')[1]) + "\t" + str(round(float(common)/float(expected), 2)) + "\n")
    else:
        w.write(con.split("/")[-1].rstrip('.' + con.split('.')[1]) + "\t0.0\n")
else:
    w.write(con.split("/")[-1].rstrip('.' + con.split('.')[1]) + "\t" + str(common) + "//" + str(expected) + "\n")
w.close()
print (str(common) + "//" + str(expected))