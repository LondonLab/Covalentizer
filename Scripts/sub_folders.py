import os
dirs = [d for d in os.listdir('./') if os.path.isdir(d)]
for d in dirs:
    print d
    os.chdir(d)
    with open('filtered_res.txt', 'r') as f:
        lines = [line.split() for line in f]
    cysteines = set([(line[1], line[2]) for line in lines])
    for i,cys in enumerate(cysteines):
        os.mkdir("CYS" + str(i))
        os.chdir("CYS" + str(i))
        with open('res.txt', 'w') as f:
            for line in lines:
                if line[1] == cys[0] and line[2] == cys[1]:
                    f.write('\t'.join(line) + '\n')
        os.chdir('../')
    os.chdir('../')
