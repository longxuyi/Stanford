import sys
import time


# print usage message if missing arguments
if __name__ == "__main__":
    if len(sys.argv) <= 3:
        print("Usage:")
        print(" $ python3 similarity.py <data_file> <output_file> <user_thresh>")
        sys.exit(0)

#create function to generate a movie dictionary from the movie data file
def create_dict(filename):
    
    #read read file and create a list of reads
    f = open(filename,'r')
    users = set()

    #create a dictionary for all movies
    dictionary = dict()

    for line in f:
        u, m, r= int(line.split()[0]), int(line.split()[1]), int(line.split()[2])
    
        if m in dictionary:
    
            dictionary[m][u] = r

        else: 
            d = {u:r}
            dictionary[m] = d
           
    f.close()
    return dictionary

#create function to count lines
def lineCount(filename):
    with open(filename,'r') as f:
        num_lines = len(f.readlines())
    
    return num_lines

#create function to count number of users
def userCount(filename):
    #read read file and create a list of reads
    f = open(filename,'r')
    users = set()
    for line in f:
        u = int(line.split()[0])
        if u not in users:
            users.add(u)
    f.close()

    return len(users)

#create function to find common users between two movies
def comUsers(movie1, movie2):

    user_set_i = set(movies[movie1].keys())
    user_set_j = set(movies[movie2].keys())

    common_users = user_set_i.intersection(user_set_j)

    return common_users    


movies = create_dict(sys.argv[1])

num_lines = lineCount(sys.argv[1])

num_users = userCount(sys.argv[1])

#convert user inputs to strings
input_file_name = str((sys.argv[1]))
onput_file_name = str((sys.argv[2]))
user_thresh = int((sys.argv[3]))

#create a list of movie ids in ascedning order
movie_id = sorted(list(movies.keys()))

# create a output file 
f = open(sys.argv[2], 'w')

#Start the timer
start = time.time()

#find the most similar movie for each movie and calculate similarity coefficient
for i in range(len(movies)):

    #initialize variables for similarity calculation
    p_max = -2
    most_similar_movie = 0
    com_users = 0
    
    #compare the movie to other movies
    for j in range(len(movies)):

        #avoid comparing movie to iteself
        if j != i: 
            
            #find a set of common users between two movies
            common_users = comUsers(movie_id[i], movie_id[j])
            
            #find number of common users
            num_comUsers = len(common_users)

        
            #if there are enough common users, calculate similarity coeeficient
            if (num_comUsers >= user_thresh):

                #initiate values for for movie a and b
                a = 0
                b = 0

                #find ratings for two movies from every common user
                for user in common_users:

                    r_ai = movies[movie_id[i]][user]
                    r_bi = movies[movie_id[j]][user]

                    a += r_ai
                    b += r_bi   


                #calculate avg ratings for movie a and b
                avg_a = a/len(common_users)
                avg_b = b/len(common_users)

                
                #calculate numerator and denominator
                numer = 0
                denom_p1 = 0
                denom_p2 = 0

                for user in common_users:
                    r_ai = movies[movie_id[i]][user]
                    r_bi = movies[movie_id[j]][user]
                    
                    numer += (r_ai - avg_a)*(r_bi - avg_b)
                    denom_p1 += (r_ai - avg_a)**2
                    denom_p2 += (r_bi - avg_b)**2
                
                denom = (denom_p1*denom_p2)**0.5

                if denom != 0:
                    p = round(numer/denom, 2)
                

                #when denominator is 0, numerator is also 0, thus no corrolation between two movies   
                else:
                    p = 0
                
                #update parameters for the most similar movie
                if p >= p_max:
                    p_max = p
                
                    most_similar_movie = movie_id[j]

                    com_users = num_comUsers

    if com_users >= user_thresh:
        f.write('{} ({},{},{})\n'.format(movie_id[i], most_similar_movie,p_max,com_users))
    else:
        f.write('{}\n'.format(movie_id[i]))

    
f.close()

#Calculate the eplapsed time
time = time.time() - start
time = round(time, 3)

#Program output in terminal
print('Input MovieLens file: {}'.format(input_file_name))
print('Output MovieLens file: {}'.format(onput_file_name))
print("Minimum Number of common users: {}".format(user_thresh))
print("Read {} lines with total of {} movies and {} users".format(num_lines, len(movies), num_users))
print("Computed similarities in {} seconds".format(time))
