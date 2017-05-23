import tensorflow as tf
import numpy as np

num_hashes = 400
num_buckets = 2**32
cpu = "/cpu:0"


def create_sparse_vec(word_list):
    num_words = len(word_list)
    indices = [[xi, 0, yi] for xi,x in enumerate(word_list) for yi,y in enumerate(x)]
    chars = list(''.join(word_list))
    return tf.SparseTensorValue(indices, chars, [num_words,1,1])



minhash = range(num_hashes)

with tf.device(cpu):
    kmers = tf.placeholder( tf.string, shape=[None], name="kmers")
    copy_kmers = kmers
    original_hashed_kmers = tf.string_to_hash_bucket_fast(kmers, num_buckets, name="hashed_kmers")
    hashed_kmers = original_hashed_kmers
    for i in range(num_hashes-1):
        minhash[i] = tf.argmin( hashed_kmers )
        hashed_kmers = tf.scalar_mul( 0x1234567887654321, hashed_kmers)
        hashed_kmers = 0x1234567887654321 + hashed_kmers
    minhash[num_hashes-1] = tf.argmin( hashed_kmers )
 
def MinHash(A):
    with tf.Session(config=tf.ConfigProto(log_device_placement=False)) as sess:
        my_minhash, my_original_hashed_kmers  = sess.run( [minhash, original_hashed_kmers],  feed_dict={kmers: A} )
    return my_original_hashed_kmers[my_minhash]

def MinHasht(A):
    with tf.Session(config=tf.ConfigProto(log_device_placement=False)) as sess:
        my_minhash,copykmers = sess.run( [minhash,copy_kmers],  feed_dict={kmers: A} )
    return copykmers[my_minhash]


def Jaccard(A,B):
    hA = MinHash(A)
    hB = MinHash(B)
    return sum(hA==hB)*1.0/num_hashes

def Jaccardt(A,B,t):
    s = 0
    hA = MinHasht(A.keys())
    hB = MinHasht(B.keys())
    for i in range(num_hashes) :
        if (hA[i] == hB[i]) and ( abs(A[hA[i]] - B[hB[i]]) < t ):
            s += 1
    return s*1.0/num_hashes     

def CreateWindow( wl ):
    window = ""
    r = np.random.randint(4, size=wl)
    for i in range(wl):
        if r[i] == 0:
            window += 'G'
        elif r[i] == 1:
            window += 'T'
        elif r[i] == 2:
            window += 'C'
        else:
            window += 'A' 
    return window

def CreateKMERPairs ( w, k ):
    kmerpairs = {}
    for i in range ( len(w) - k):
        for j in range (i+1,len(w) - k ):
            a = w[i:i+k]
            b = w[j:j+k]
            d = j - i
            ab= a+b
            if ab in kmerpairs:
                kmerpairs[ab] = min(d,kmerpairs[ab])
            else :
                kmerpairs[ab] = d
    return kmerpairs

def CreateKMERPairsThreshold ( kmers, t ):
    newkmers = {}
    for k in kmers.keys():
        if kmers[k] < t:
            newkmers[k+"0"] = 0
        else :
            newkmers[k+"1"] = 1
    return newkmers




w1 = CreateWindow( 500 )
kmers1 = CreateKMERPairs ( w1, 8 )
w2 = CreateWindow( 500 )
kmers2 = CreateKMERPairs ( w2, 8 )
kmers1t = CreateKMERPairsThreshold ( kmers1, 5 )
kmers2t = CreateKMERPairsThreshold ( kmers2, 5 )
print Jaccard(kmers1t.keys(),kmers2t.keys())
print Jaccardt(kmers1, kmers2, 10 )