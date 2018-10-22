in_f = open ('ERK_data.txt')
out_f = open ('ERK_filtered_data.txt', 'w')
for line in in_f:
    l = line.split ()
    if '#' in l:
        end_idx = l.index ('#')
    else:
        end_idx = len (l)

    for i in range (end_idx):
        if i > 0 and l[i].replace ('.', '', 1).isdecimal ():
            l[i] = str (round(float (l[i]) * 100, 2))

    out_f.write (' '.join (l) + '\n')


