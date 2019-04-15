import os
import csv
import numpy as np
import operator
def parseVCF(vcf, path, filter = 5):
    vcf_array = []
    GT_calls = {'0/0': '0', '1/1': '1', '0/1': '1', './.': '0'}
    with open(os.path.join(path, vcf), 'r') as myVCF:
        vcfReader = csv.reader(myVCF, delimiter='\t')
        for field in vcfReader:
            if field[0].startswith('#'):
                pass
            else:
                TE_info = field[7] #Info on the TE
                ac = [True if filt == 'ac0' else False for filt in field[6].split(';')]

               # print(TE_info)
                if int(TE_info.split(';')[1].split('=')[1]) >= filter and any(ac) == False: #filter by the assessment score
                    chrom = field[0]
                    indivs = field[9:24]
                    position = field[1]
                    GTs = []
                    for gt_call in indivs:
                        mei = GT_calls[gt_call.split(':')[0]]
                        GTs.append(mei)
                    mei_line = [chrom, position] + GTs
                    vcf_array.append(mei_line)
                else:
                    pass
        myVCF.close()

        return vcf_array
def write_out(path, vcf, vcf_array):
    ##write out GT calls:
    outPUT = os.path.join(path, vcf.split('.')[0]+'mei.tsv')
    strains = ['BL6.NTac', 'BL6.NCrl', 'BL6.NHsd', 'BL6N-TyrCBrdCrlCrl', 'BL6.JBomTac', 'BL10.ScSnJ', 'BL6.J', 'BL10.SnJ', 'BL10.ScCr', 'BL6.NJ', 'BL10.J', 'BL6.ByJ', 'BL6.JEiJ', 'BL10.ScNHsd']
    tip_order = ["BL6.NTac", "BL6.NJ", "BL6.NCrl", "BL6N-TyrCBrdCrlCrl", "BL6.NHsd", "BL6.ByJ", "BL6.JBomTac", "BL6.J", "BL6.JEiJ", "BL10.ScSnJ", "BL10.SnJ", "BL10.ScCr", "BL10.ScNHsd", "BL10.J"]


    new_order = {} #getting some weird errors
    for tip in range(len(tip_order)):
        new_order[tip_order[tip]] = tip


    with open(outPUT, 'w') as myTSV:

        header = "chr\tpos"
        for mouse in tip_order:
            header += '\t'+mouse
        header = header + '\n'
        myTSV.write(header)

        for mei in vcf_array:
            delim = '\t'
            mei_vector = [np.nan for x in range(14)]
            mei_zip = zip(strains, mei[2:17])
            for item in mei_zip:
                mei_vector[new_order[item[0]]] = item[1]   #create new list with the new indices
            newLine = delim.join(mei[0:2] + mei_vector) + '\n'
            myTSV.write(newLine)
        myTSV.close()
def write_trans(path, vcf, vcf_array):
    outPUT = os.path.join(path, vcf.split('.')[0] + 'mei.tsv')
    strains = ['BL6.NTac', 'BL6.NCrl', 'BL6.NHsd', 'BL6N-TyrCBrdCrlCrl', 'BL6.JBomTac', 'BL10.ScSnJ', 'BL6.J', 'BL10.SnJ', 'BL10.ScCr', 'BL6.NJ', 'BL10.J', 'BL6.ByJ', 'BL6.JEiJ', 'BL10.ScNHsd']
    tip_order = ["BL6.NTac", "BL6.NJ", "BL6.NCrl", "BL6N-TyrCBrdCrlCrl", "BL6.NHsd", "BL6.ByJ", "BL6.JBomTac", "BL6.J", "BL6.JEiJ", "BL10.ScSnJ", "BL10.SnJ", "BL10.ScCr", "BL10.ScNHsd", "BL10.J"]
    new_order = {}  # getting some weird errors
    for tip in range(14):
        new_order[strains[tip]] = tip_order.index(strains[tip])


    with open(outPUT, 'w') as myOut:
        all_cols = ''
        for posHeader in vcf_array:
            position = '_'.join(posHeader[0:2])
            all_cols += '\t' + position
        myOut.write('strain{0}\n'.format(all_cols))



        #write GT calls
        correct_order = sorted(new_order.items(), key=operator.itemgetter(1))
        for strain in correct_order:
            strain_ind = strains.index(strain[0]) + 2

            strain_gt = ''
            for gt in vcf_array:
                strain_gt += '\t' + gt[strain_ind]
            myOut.write('{0}{1}\n'.format(strain[0], strain_gt))
        myOut.close()


for file in os.listdir('/home/iskander/Documents/Clark_lab/VCF_outputs'):
    if file.endswith('.vcf'):
        arr = parseVCF(path='/home/iskander/Documents/Clark_lab/VCF_outputs', vcf=file)

        write_trans(path='/home/iskander/Documents/Clark_lab/VCF_outputs', vcf=file, vcf_array=arr)



