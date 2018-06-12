
# coding: utf-8

# In[6]:


##Yousuf's Method. The input is meant to be a list of equal-sized k-mers and the goal is figure out the total of number of 
## nucleotide counts at each posiiton and return that as a dictionary of lists (each list is a different nucleotide)
def Count(Motifs):
    #I am creating a dictionary of nucleotides which will contain a list of nt frequency at each position
    count = {"A": [], "C": [], "G": [], "T": []}
    #d represents the number of different sequences I want to find the motif in
    d = len(Motifs)
    #k represents the length of the kmer I'm looking at... since all of them are the same length I'm just looking at
    #the first one
    k = len(Motifs[0])
    #I am first iterating through the positions of each motif, starting at position zero and ending with the final position k
    for i in range(k):
        A = 0
        T = 0
        C = 0
        G = 0
        #I am iterating through every string at position i and will add up how many nucleotides go to what then append the 
        #number to the list
        for kmer_number,kmer in enumerate(Motifs):
            if kmer[i] == "A":
                A +=1
            if kmer[i] == "T":
                T +=1
            if kmer[i] == "C":
                C +=1
            if kmer[i] == "G":
                G+=1
            if kmer_number == len(Motifs)-1:
                count["A"].append(A)
                count["T"].append(T)
                count["C"].append(C)
                count["G"].append(G)
    return count


# In[7]:


##Using the hints/advice given to me by the program, I will design this second 

def Count2(Motifs):
    #Making an empty dictionary
    count = {}
    #Creating the dictionary with the appropriate number of zeros in every position for every nucleotide already inputted
    #k is th length of the motif. Each motif is the same length so we can just use the first one
    #The list iterates through the 4 nucleotides of life, ACGT and creates an empty list [] with the key being the symbol
    #Once an empty list is created with the key as the nucleotide symbol, a list of zeros is added
    #The number of zeros corresponds to k
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
             count[symbol].append(0)
    #t is the number of motifs we need t count
    t = len(Motifs)
    #iterating from 0 the the max number of sequences
    for i in range(t):
        #For each motif, we will iterate through the length
        for j in range(k):
            #We are going to read what the nucleotide is at the motif we're looking at --> i and what position --> j
            symbol = Motifs[i][j]
            #add +1 to the dictionary 
            count[symbol][j] += 1
    return count


# In[8]:


#This is a program that will normalize the nucleotide counts at each position to a decimal frequency 
def Profile(Motifs):
    #number of motifs we have to count
    t = len(Motifs)
    #length of each motif..all of them are the same length so we just count the length of the first one
    k = len(Motifs[0])
    #getting our raw nucleotide counts from the count program
    count = Count(Motifs)
    #setting up a new dictionary using the technique from count2 where all the positions are 0's
    normalized_count = {}
    for symbol in "ACGT":
        normalized_count[symbol] = []
        for j in range(k):
             normalized_count[symbol].append(0)
    #iterating through all the positions starting from 0
    for position in range(k):
        total = 0
        #count the total number of nucleotide counts at each position
        for nt in "ACGT":
            total += count[nt][position]
        #divide the raw nucleotide total by the total number of all nucleotides
        for nt in "ACGT":
            normalized_count[nt][position] = float(count[nt][position])/total
    return normalized_count
        


# In[9]:


#This is a program that will give me back a consensus sequence based on count
def Consensus(Motifs):
    consensus = ""
    k = len(Motifs[0])
    count = Count(Motifs)
    for position in range(k):
        max = 0
        max_nt = ''
        for nt in "ACGT":
            if count[nt][position] > max:
                max = count[nt][position]
                max_nt = nt
        consensus += max_nt
    return consensus


# In[10]:


#This is a program that will score the consensus sequence based on how many nucleotides DO NOT match the consensus sequence
def Score(Motifs):
    consensus = Consensus(Motifs)
    count = Count(Motifs)
    k = len(Motifs[0])
    score = 0
    for position in range(k):
        current_nt = consensus[position]
        for symbol in "ACGT":
            if symbol != current_nt:
                score += count[symbol][position]
    return score
        


# In[11]:


#This program will calculate the probability of obtaining a k-mer given the output(probability matrix) from Profile 
#See Biology Meets Programming - Motif Counter if confused
def Pr(Text, Profile):
    k = len(Text)
    probability = 1
    matrix = Profile
    for position in range(k):
        probability*= Profile[Text[position]][position]
    return probability


# In[12]:


#This program will find the MOST likely k-mer given the: k-mer size, the profile matrix, and sequence we want to search through
def ProfileMostProbablePattern(text, k, profile):
    #Since this sequence is circular, we need to add a bit of the beginning of the sequence to end for long motifs
    extendedGenome = text + text[0:len(text)//2]
    #just initializing some variables
    best_prob = -1
    best_kmer = extendedGenome[0:k]
    #So first I'll be iterating through each position
    for i in range(len(text)-k+1):
        #Checks probability of k-mer i to i+k
        current_prob = Pr(extendedGenome[i:i+k], profile)
        #if probability is better than previous one, replace.
        if current_prob > best_prob:
            best_prob = current_prob
            best_kmer = extendedGenome[i:i+k]
    return best_kmer


# In[13]:


#This program will find the best k-mer given: list of sequence strings, size of k-mer, number of sequences to search through 
def GreedyMotifSearch(Dna, k, t):
    #The best motif is set as the first k-mer of each sequence, so you have a list of t-number of k-mers
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    #Now n is the length of the first sequence I have
    n = len(Dna[0])
    #i will be the starting position of where I will be drawing my k-mers from
    for i in range(n-k+1):
        #Empty list creation
        Motifs = []
        #I add a k-mer of appropriate length from the first sequence string to the list Motifs
        Motifs.append(Dna[0][i:i+k])
        #iterating from 0 to the number of sequences you have
        for j in range(1, t):
            #P is equal to the probability matrix generated by profile, which will count how similar nucleotides are and return
            #it as a percentage of A's or T's at a position. This will only have up to t-number of k-mers in it
            P = Profile(Motifs[0:j])
            #Add the most similar k-mer from the following sequences 1,2,3...t to this motifs list
            Motifs.append(ProfileMostProbablePattern(Dna[j], k, P))
            #Now we will check to see if the score of our generated 5 k-mers is better than the previous one. Thus we end up
            #eventually with the 5 best k-mers 
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs
    


# In[14]:





# In[15]:


####The code below this is for revising the previous greedy motif method into a better one


# In[16]:


#This is a revised method of counting, since the original method suffers from Cromwell's problem i.e. assigning something a 
#probability of zero is not the best way to create a set... everything should have a nonzero percent chance of happening,
#no matter how unlikely
def CountWithPseudocounts(Motifs):
    #Making an empty dictionary
    count = {}
    #Creating the dictionary with the appropriate number of zeros in every position for every nucleotide already inputted
    #k is th length of the motif. Each motif is the same length so we can just use the first one
    #The list iterates through the 4 nucleotides of life, ACGT and creates an empty list [] with the key being the symbol
    #Once an empty list is created with the key as the nucleotide symbol, a list of zeros is added
    #The number of zeros corresponds to k
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
             count[symbol].append(1)
    #t is the number of motifs we need t count
    t = len(Motifs)
    #iterating from 0 the the max number of sequences
    for i in range(t):
        #For each motif, we will iterate through the length
        for j in range(k):
            #We are going to read what the nucleotide is at the motif we're looking at --> i and what position --> j
            symbol = Motifs[i][j]
            #add +1 to the dictionary 
            count[symbol][j] += 1
    return count    


# In[17]:


#Revised ProfileCount
def ProfileWithPseudocounts(Motifs):
    #number of motifs we have to count
    t = len(Motifs)
    #length of each motif..all of them are the same length so we just count the length of the first one
    k = len(Motifs[0])
    #getting our raw nucleotide counts from the count program
    count = CountWithPseudocounts(Motifs)
    #setting up a new dictionary using the technique from count2 where all the positions are 0's
    normalized_count = {}
    for symbol in "ACGT":
        normalized_count[symbol] = []
        for j in range(k):
             normalized_count[symbol].append(0)
    #iterating through all the positions starting from 0
    for position in range(k):
        total = 0
        #count the total number of nucleotide counts at each position
        for nt in "ACGT":
            total += count[nt][position]
        #divide the raw nucleotide total by the total number of all nucleotides
        for nt in "ACGT":
            normalized_count[nt][position] = float(count[nt][position])/total
    return normalized_count
     


# In[18]:


#Updated GreedyCount 
def GreedyMotifSearchWithPseudocounts(Dna, k, t):

    #The best motif is set as the first k-mer of each sequence, so you have a list of t-number of k-mers
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    #Now n is the length of the first sequence I have
    n = len(Dna[0])
    #i will be the starting position of where I will be drawing my k-mers from
    for i in range(n-k+1):
        #Empty list creation
        Motifs = []
        #I add a k-mer of appropriate length from the first sequence string to the list Motifs
        Motifs.append(Dna[0][i:i+k])
        #iterating from 0 to the number of sequences you have
        for j in range(1, t):
            #P is equal to the probability matrix generated by profile, which will count how similar nucleotides are and return
            #it as a percentage of A's or T's at a position. This will only have up to t-number of k-mers in it
            P = ProfileWithPseudocounts(Motifs[0:j])
            #Add the most similar k-mer from the following sequences 1,2,3...t to this motifs list
            Motifs.append(ProfileMostProbablePattern(Dna[j], k, P))
            #Now we will check to see if the score of our generated 5 k-mers is better than the previous one. Thus we end up
            #eventually with the 5 best k-mers 
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs


# In[19]:


#This is a new way to find motifs in sequences


# In[37]:


#Find the best k-mer from each sequence of DNA present in the list DNA
def Motifs(Profile, Dna):
    kmerList = []
    for seq in Dna:
        kmerList.append(ProfileMostProbablePattern(seq, len(Profile['A']), Profile))
    return kmerList


# In[39]:


#RandomMotifs(Dna, k, t) uses 
#random.randint to choose a random k-mer from each of t different strings Dna, and returns a list of t strings.
def RandomMotifs(Dna, k, t):
    import random
    finalKmer = []
    for seq in Dna:
        tempKmer = []
        for i in range(len(seq)-k+1):
            tempKmer.append(seq[i:i+k])
        finalKmer.append(tempKmer[random.randint(0, len(tempKmer)-1)])
    return finalKmer
        
        



# In[40]:


#This function will keep finding random motifs that are better than the previous one for N-runs
def RandomizedMotifSearch(Dna, k, t, N):
    M = RandomMotifs(Dna, k, t)
    BestMotifs = M
    runs = 0
    #A list of random Motifs has been generated
    while True & runs < N:
        #Generate a profile based on the motifs
        runs +=1
        Profile = ProfileWithPseudocounts(M)
        M = Motifs(Profile, Dna)
        if Score(M) < Score(BestMotifs):
            BestMotifs = M
        else:
            return BestMotifs 


# In[ ]:


##The following functions are meant to improve the randomized counting method 


# In[50]:


#This function takes a set of probabilities and normalizes them so they add up to 1
def Normalize(Probabilities):
    sum = 0
    normProb = {}
    for key in Probabilities:
        normProb[key] = 0
        sum += Probabilities[key]
    factor = 1.0/sum
    for key in normProb:
        normProb[key] = Probabilities[key] * factor
    return normProb
    


# In[90]:


#return a random k-mer given a dictionary of k-mers and probabilities
def WeightedDie(Probabilities):
    import random
    newProb = Normalize(Probabilities)
    kmer = '' 
    dieRoll = random.uniform(0, 1)
    iterableSum = 0
    for seq in newProb:
        iterableSum += newProb[seq]
        if iterableSum >= dieRoll:
            kmer = seq
            break
    return kmer


# In[92]:


def ProfileGeneratedString(Text, profile, k):
    import random
    n = len(Text)
    probabilities = {} 
    for i in range(0,n-k+1):
        probabilities[Text[i:i+k]] = Pr(Text[i:i+k], profile)
    return WeightedDie(probabilities)


# In[94]:


#The final program
def GibbsSampler(Dna, k, t, N):
    import random
    BestMotifs = []
    Motifs = RandomMotifs(Dna, k, t)
    BestMotifs = Motifs
    for j in range(N):
        i = rand.int(1, t)
        
        Mi = Motifs[i]
        
        del Motifs[i]
        
        prof = ProfileWithPseudocounts(Motifs)
        
        Motifs = Motifs(prof,Mi)
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs
    

