#Get SHAPE data
structure = stringToList(structure)
structure = numpy.array(structure)
structure = structure.astype('float64')
#Get RNA sequence
sequence= protein_2_embed(sequence)
#Merge
seq_str = np.vstack((sequence,structure))