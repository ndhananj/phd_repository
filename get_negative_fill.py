################################################################################
# Get the negative fill of a bottlneck
# Originally made by Nithin Dhananjayan (ndhanananj@ucdavis.edu)
# Usage : python <this_file_name> <bottleneck_file> <negative_fill_file>
# example : python get_negative_fill.py bottleneck.pdb fill.pdb
################################################################################


from projections import *

if __name__ == '__main__':
    get_negative_fill_3D(sys.argv[1],sys.argv[2])
