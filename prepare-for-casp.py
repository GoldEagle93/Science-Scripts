with open('target.cgi.txt', 'r') as inputfile:
    temp = inputfile.readlines()
inputfile.close()

indices = ['0', '1', '2', '3', '4']
with open('submission.pdb', 'w') as outputfile:
    for index in indices:
        with open('unique.c%s.pdb'%index, 'r') as inputfile:
            cluster = inputfile.readlines()
        inputfile.close()

        with open('back.%s.apf'%index, 'r') as inputfile:
            bfactors = inputfile.readlines()
        inputfile.close()

        bf = {}

        for line in bfactors[1:]:
            bf[line.split()[0][:-4]] = line.split()[1][:4]

        model = []
        for i in cluster:
            if '           H  ' not in i:
                model.append(i)

        outfile = []
        for x in temp:
            xx = x.split()
        # for x in temp[8200:8230]:
            if xx[0] == 'TER':
                outfile.append(x)
                continue
            for y in model:
                yy = y.split()
                if yy[-1] == 'H' or yy[0] == 'TER' or 'END' in yy[0] or 'CRYST' in y or  'END' in xx:
                    continue
                elif xx[4] == '0' and int(xx[5])+0 == int(yy[4]) and xx[2] == yy[2]:
                    outfile.append(
                    x[:30] + y[30:56] + x[56:60] + '  ' + bf[yy[4]][:4] + x[66:]
                    )
        outputfile.write('MODEL %s\n'%str(int(index)+1))
        outputfile.write('PARENT N/A\n')
        outputfile.write('ATOM      1  P     G 0   1       0.000   0.000   0.000  1.00  0.00           P\n')
        outputfile.write('ATOM      2  OP1   G 0   1       0.000   0.000   0.000  1.00  0.00           O\n')
        outputfile.write('ATOM      3  OP2   G 0   1       0.000   0.000   0.000  1.00  0.00           O\n')
        for line in outfile:
            outputfile.write(line)
        outputfile.write('ENDMDL\n')
    outputfile.write('END\n')
outputfile.close()
