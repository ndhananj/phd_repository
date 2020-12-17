from modes import *

import sys

def embedding_comparison(dirs,labels):
    em_files = [dir+'/embedding.npy' for dir in dirs]
    fig, ax = plt.subplots()
    for i in range(len(em_files)):
        prob=load_matrix(em_files[i])
        ax.scatter([i+1 for i in range(len(prob))],prob,label=labels[i])
    plt.xlabel("Number of modes in space")
    plt.ylabel('"Overlap"')
    plt.legend()
    plt.show()

if __name__ == '__main__':
    dirs = sys.argv[1::2]
    labels = sys.argv[2::2]
    embedding_comparison(dirs,labels)
