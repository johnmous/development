import sys
import click
import os
import glob
from pathlib import Path

@click.command()
@click.option('--srr2sexweekfile', type=click.Path(exists=True, readable=True),
              required=True, help='Table with srrid (column 1), \
breakpoint for SRRs with more than one embryo-that is the cell number where the \
switch between embryos happens (column 5), embryo ID (column 6)')
@click.option('--srrid', type=str, help='SRR ID of current file')
@click.option('--barcd2cellidfile', type=click.Path(exists=True, readable=True),
required=True, help='File with two fields, srrid and Cell ID. Provided by the authors')
@click.option('--path2barcodes', type=click.Path(exists=True,
readable=True), required=True, help='Path with alignment files, file name is cell barcode')
@click.option('--path2cellids', type=click.Path(writable=True),
required=True, help='Path with alignment files, file name is cell ID')
def main(barcd2cellidfile, path2barcodes, path2cellids, srr2sexweekfile, srrid):
    SRR2brkpnt2sexWeekArray = SRR(srr2sexweekfile)
    barcd2CellId = open(barcd2cellidfile, "r")
    file_and_path = glob.glob(path2barcodes+"/*.bam")
    experimentBarcd =  [os.path.basename(f) for f in file_and_path]
    if not os.path.exists(path2cellids):
        os.mkdir(path2cellids)

    barcd2sc ={}
    for line in barcd2CellId:
        l = line.split()
        barcode = l[0]
        sc = l[1]
        scInt = int(sc[2:])
        barcd2sc[barcode] = [scInt, sc]

    for barcode_bam in experimentBarcd:
        barcode = barcode_bam[:-4]
        breakpoint = SRR2brkpnt2sexWeekArray[srrid][0]
        sexWeekArray = SRR2brkpnt2sexWeekArray[srrid][1]
        scInt = barcd2sc[barcode][0]
        scStr = barcd2sc[barcode][1]
        if scInt >= breakpoint:
            sample = sexWeekArray[-1]
        else:
            sample = sexWeekArray[0]
        os.symlink(path2barcodes+'/'+barcode+'.bam',
                    path2cellids+'/'+sample+"_"+scStr+'.bam')
    barcd2CellId.close()
    Path(path2cellids+'/'+srrid+'_createlinks_done.txt').touch()


def SRR(srr2sexweekfile):
    SRR2sexWeek = open(srr2sexweekfile, 'r')
    next(SRR2sexWeek)
    SRR2brkpnt2sexWeekArray = {}
    for line in SRR2sexWeek:
        splitLine = line.split()
        if len(splitLine) != 6:
            sys.exit("Table should have exactly 6 columns")
        srrid = splitLine[0]
        breakpoint = int(splitLine[4])
        if breakpoint != 0:
            sexWeekArray = splitLine[5].split(",")
        else:
            sexWeekArray = [splitLine[5]]
        SRR2brkpnt2sexWeekArray[srrid] = [breakpoint, sexWeekArray]
    SRR2sexWeek.close()
    return(SRR2brkpnt2sexWeekArray)

if __name__ == '__main__':
    main()
