f_m = open("/Volumes/ATUL_6TB/Work/Projects/Microglia_DEG/SNP_ID.txt", 'w', 1)
with open ("/Volumes/ATUL_6TB/Work/Projects/Microglia_DEG/6m_FADKOvsFAD_Ortholog.txt", 'r') as myFile:
    line = myFile.readline ()
    for line in myFile:
        line_list = line.split("\t")
        l = []
        s = int (line_list [12]) - 1000000
        e = int (line_list [13]) + 1000001
        for i in range (s, e):
            SNP = str (line_list [11]) + ":" + str(i)
            l.append (SNP)
        for id in l:
            f_m.write (id + "\n")
