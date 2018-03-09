import os
command = 'ls *.pdf > tmp.dat'
os.system(command)
f = open('tmp.dat')
for line in f:
    line = line.split('.pdf')
    command = 'cp %s.pdf eps2:%s.eps '%(line[0],line[0])
    os.system(command)
    #print command

