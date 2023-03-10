'''Some simple tools to handle single cell matrices in python.
Just more comfortable for me than R, definitely slower and worse
in most senses.'''

import os
import subprocess
import pandas as pd

'''outputs a cell-indexed (rows) into a sparse pandas dataframe of matrix.'''
def MatrixToPanda(cellBarcodes, geneNames, matrix, cellHeader, geneHeader):
    #read-in cells
    col = ['gene', 'cell', 'UMI']
    output = []
    #no header in non-annotated file
    cells = []
    with open(cellBarcodes) as infile:
        if cellHeader:
            infile.readline()
        for line in infile:
            cells.append(line[:-1].split('\t')[0])

    #read-in genes.
    genes = []
    with open(geneNames) as infile:
        if geneHeader:
            infile.readline()
        for line in infile:
            genes.append(line[:-1].split('\t')[0])
    #now read in matrix, need to ignore first few lines
    with open(matrix) as infile:
        skiplines = 3
        while skiplines > 0:
            infile.readline()
            skiplines -= 1
        for line in infile:
            values = line[:-1].split(' ')
            values[0] = genes[int(values[0]) - 1]
            values[1] = cells[int(values[1]) - 1]
            values[2] = int(values[2])
            output.append(values)
    output = pd.DataFrame(output, columns = col)
    return output

'''outputs a sparse matrix file and associated cell and barcode files from a tidy pandas matrix. Keep original gene list to keep 0-value genes.
Organize genes and cells by alphabetical order to speed things. Depreciated for now...takes too long'''
def PandaToMatrix(indata, geneFile, outCells, outMatrix, geneHeader):
    #read in genes
    geneList = []
    with open(geneFile) as infile:
        if geneHeader:
            infile.readline()
        for line in infile:
            geneList.append(line[:-1].split('\t')[0])
    cellList = indata.cell.unique().tolist()
    #now change tidy dataframe into sparse matrix

    print('changing cells')
    indata.cells = indata.apply(lambda row: cellList.index(row.cell) + 1, axis = 1)
    print('changing genes')
    indata.genes = indata.apply(lambda row: geneList.index(row.gene) + 1, axis = 1)
    with open(outMatrix, 'w') as matrixFile, open(outCells, 'w') as cellFile:
        for cell in cellList:
            cellFile.write(cell + '\n')
        #first, write header for matrix file
        matrixFile.write('%%MatrixMarket matrix coordinate integer general\n%\n')
        #for a matrix, first column is genes, second cells, and third counts. Top of matrix is summation of different genes, different cells, and total counts
        geneNum = str(len(geneList))
        cellNum = str(len(cellList))
        UMIcount = str(indata.UMI.sum())
        matrixFile.write(" ".join([geneNum, cellNum, UMIcount]) + '\n')
        matrixFile.write(indata.to_csv(sep =' ', header=False, index = False))


def storeMatrixPanda(inMatrix):
    col = ['gene', 'cell', 'UMI']
    output = []

    with open(inMatrix) as infile:
        skiplines = 3
        while skiplines > 0:
            infile.readline()
            skiplines -= 1
        for line in infile:
            values = line[:-1].split(' ')
            values[0] = int(values[0])
            values[1] = int(values[1])
            values[2] = int(values[2])
            output.append(values)
    output = pd.DataFrame(output, columns = col)
    return output

#outputs new matrix. Assumes only editing out cells, not genes. Iterate through cells ONCE, generating a new linkage with old numbers. Then lookup for number
#changing becomes trivial
def sparsePandaToMatrix(originalCells, inPandas, outMatrix, outCells, geneHeader, geneFile):

    #no header in non-annotated file
    cells = []
    #enforce integer values
    inPandas['cell'] = inPandas.cell.astype(int)
    inPandas['gene'] = inPandas.gene.astype(int)
    inPandas['UMI'] = inPandas.UMI.astype(int)
    with open(originalCells) as infile:
        for line in infile:
            cells.append(line[:-1].split('\t')[0])
    #this will be the new order
    newCells = inPandas.cell.unique().tolist()
    cellsToWrite = []
    #new order of cells
    for cell in newCells:
        cellsToWrite.append(cells[cell-1])
    #modify old list to reference new indices, will work on a 0 indexed system and increment in output
    for position, cell in enumerate(cells):
        try:
            cells[position] = cellsToWrite.index(cell)
        except:
            cells[position] = 0
    inPandas.cell = inPandas.apply(lambda row: cells[int(row.cell - 1)] + 1, axis = 1)
    with open(outMatrix, 'w') as matrixFile, open(outCells, 'w') as cellFile:
        for cell in cellsToWrite:
            cellFile.write(cell + '\n')
        #first, write header for matrix file
        matrixFile.write('%%MatrixMarket matrix coordinate integer general\n%\n')
        #for a matrix, first column is genes, second cells, and third counts. Top of matrix is summation of different genes, different cells, and total counts
        #get number of genes. Easy to do, number of lines in genefile unless geneHeader, then that number -1
        proc = subprocess.Popen(" ".join(['wc -l', geneFile]), shell=True, stdout=subprocess.PIPE)
        geneNum = int(proc.stdout.readline().decode("utf-8").split(' ')[0])
        if geneHeader:
            geneNum -= 1
        cellNum = str(len(newCells))
        UMIcount = str(int(inPandas.UMI.sum()))
        geneNum = str(geneNum)
        matrixFile.write(" ".join([geneNum, cellNum, UMIcount]) + '\n')
        columns = ['gene', 'cell', 'UMI']
        inPandas = inPandas[columns]
        matrixFile.write(inPandas.to_csv(sep =' ', header=False, index = False))


#as above, saves time by not needing to fully annotate. Returns row-number keyed dictionary of genes.
def sparseGeneList(geneList, geneFile, geneHeader):
    geneReturn = []
    geneOrder = []
    with open(geneFile) as infile:
        if geneHeader:
            infile.readline()
        for line in infile:
            geneOrder.append(line.split('\t')[0])
    for gene in geneList:
        try:
            geneReturn.append(geneOrder.index(gene) + 1)
        except:
            pass
    return geneReturn


def MatrixToPandaSoupX(cellBarcodes, geneNames, matrix, cellHeader, geneHeader):
    #read-in cells
    col = ['gene', 'cell', 'UMI']
    output = []
    #no header in non-annotated file
    cells = []
    with open(cellBarcodes) as infile:
        if cellHeader:
            infile.readline()
        for line in infile:
            cells.append(line[:-1].split('\t')[0])

    #read-in genes.
    genes = []
    with open(geneNames) as infile:
        if geneHeader:
            infile.readline()
        for line in infile:
            genes.append(line[:-1].split('\t')[0])
    #now read in matrix, need to ignore first few lines
    with open(matrix) as infile:
        skiplines = 3
        while skiplines > 0:
            infile.readline()
            skiplines -= 1
        for line in infile:
            values = line[:-1].split(' ')
            values[0] = genes[int(values[0]) - 1]
            values[1] = cells[int(values[1]) - 1]
            values[2] = float(values[2])
            output.append(values)
    output = pd.DataFrame(output, columns = col)
    return output
