model='turtle_real/';
path=['../Models/' model];

T1 = read_curvepoints([path 'curvepoints1.txt']);
T2 = read_curvepoints([path 'curvepoints2.txt']);
T3 = read_curvepoints([path 'curvepoints3.txt']);

diag1 = norm(max(T1)-min(T1))
diag2 = norm(max(T2)-min(T2))
diag3 = norm(max(T3)-min(T3))

aver1 = 1.952498/diag1
aver2 = 2.074066/diag2
aver3 = 5.005622/diag3

max1 = 14.422205/diag1
max2 = 13.601471/diag2
max3 = 50.635956/diag3