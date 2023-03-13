import numpy as np
import pandas as pd

def load(filename, include_demographics=False):
    df = pd.read_csv(f"{filename}")
    if not include_demographics:
        df = df.drop(columns=["Demographic"])
    
    return df

def get_p_x_given_y(x_column, y_column, df):
    """
    Computes P(X = 1 | Y = 1) and P(X = 1 | Y = 0), where X is a single feature (column).
    x_column: name of the column containing the feature X.
    y_column: name of the class containing the class label.

    return: [P(X = 1 | Y = 1), P(X = 1 | Y = 0)]
    """

    ### YOUR CODE HERE
    # Hint: you can check which rows of a feature (e.g. "x2") are equal to some value (e.g. 1) by doing
    # df["x2"] == 1.
    # Hint: remember to use Laplace smoothing.


    #df[x_column].value_counts()[1]

    num_y0 =0
    num_y1 =0
    num_x1y1 = 0
    num_x1y0 = 0

    for i in range(df.shape[0]):
        if df.at[i,y_column] == 1:
            num_y1 +=1
            if df.at[i,x_column] == 1:
                num_x1y1 +=1
        if df.at[i,y_column] == 0:
            num_y0 +=1
            if df.at[i,x_column] == 1:
                num_x1y0 +=1
    #print("x1y0: {}, x1y1: {}, y0: {}, y1: {}".format(num_x1y0, num_x1y1, num_y0, num_y1))

    p_y0 = (num_x1y0+1)/(num_y0+2)
    p_y1 = (num_x1y1+1)/(num_y1+2)
    ### END OF YOUR CODE

    return [p_y0, p_y1]


def get_all_p_x_given_y(y_column, df):
    # We want to store P(X_i=1 | Y=y) in p_x_given_y[i][y]
    all_p_x_given_y = np.zeros((df.shape[1]-1, 2))

    ### YOUR CODE HERE
    # Hint: df.columns gives a list of all the columns of the DataFrame.
    # Hint: remember to skip the "Label" column.

    for i in range(all_p_x_given_y.shape[0]):
        p_y0,p_y1 = get_p_x_given_y(df.columns[i], y_column, df)
        all_p_x_given_y[i][0] = p_y0
        all_p_x_given_y[i][1] = p_y1
    ### END OF YOUR CODE

    #print(all_p_x_given_y)

    return all_p_x_given_y

def get_p_y(y_column, df):
    """
    Compute P(Y = 1)
    """
    count = 0
    for i in range(df.shape[0]):
        if df.at[i,y_column] == 1:
            count +=1

    p_y = count/df.shape[0]
    ### END OF YOUR CODE
    return p_y

    
def joint_prob(xs, y, all_p_x_given_y, p_y):
    """
    Computes the joint probability of a single row and y
    """

    ### YOUR CODE HERE
    # Hint: P(X, Y) = P(Y) * P(X_1 | Y) * P(X_2 | Y) * ... * P(X_n | Y)

    prob = 1

    if y == 1:
        prob *= p_y
    if y == 0:
        prob *= (1-p_y)

    for i in range(all_p_x_given_y.shape[0]):
        if xs[i] == 1:
            prob *= all_p_x_given_y[i][y]
        else:
            prob *= (1 - all_p_x_given_y[i][y])

    ### END OF YOUR CODE
    
    return prob


def get_prob_y_given_x(y, xs, all_p_x_given_y, p_y):
    """
    Computes the probability of a y given x.
    """

    n, _ = all_p_x_given_y.shape # n is the number of features/columns


    ### YOUR CODE HERE
    # Hint: use the joint probability function.
    prob_y = joint_prob(xs, y, all_p_x_given_y, p_y)
    prob_y_inv = joint_prob(xs, 1-y, all_p_x_given_y, p_y)
    prob_y_given_x = prob_y/(prob_y+prob_y_inv)
    
    ### END OF YOUR CODE
    return prob_y_given_x
#============================================

def compute_accuracy(all_p_x_given_y, p_y, df):
    # split the test set into X and y. The predictions should not be able to refer to the test y's.
    X_test = df.drop(columns="Label")
    y_test = df["Label"]
    

    num_correct = 0
    total = len(y_test)
    cur_class = 0

    ### YOUR CODE HERE
    # Hint: we predict 1 if P(Y=1|X) >= 0.5.
    # Hint: to loop over the rows of X_test, use:
    #       for i, xs in X_test.iterrows():
    for i, xs in X_test.iterrows():

        cur_prob = get_prob_y_given_x(1, xs, all_p_x_given_y, p_y)
        if cur_prob >= 0.5:
            cur_class = 1
        else: 
            cur_class = 0

        if cur_class == y_test[i]:
            num_correct = num_correct+1   

    accuracy = num_correct/total

    ### END OF YOUR CODE

    return accuracy

def py1_given_x(x_column, y_column, df):
    """
    compute P(Y=1|Xi =0) and P(Y=1|Xi =1)
    """
    num_x0 =0
    num_x1 =0
    num_y1x1 = 0
    num_y1x0 = 0

    for i in range(df.shape[0]):
        if df.at[i,x_column] == 1:
            num_x1 +=1
            if df.at[i,y_column] == 1:
                num_y1x1 +=1
        if df.at[i,x_column] == 0:
            num_x0 +=1
            if df.at[i,y_column] == 1:
                num_y1x0 +=1
   
    p_y1x0 = (num_y1x0+1)/(num_x0+2)
    p_y1x1 = (num_y1x1+1)/(num_x1+2)

    ratio = p_y1x1/p_y1x0

    return ratio

def three_movies(df):

    ratio_list = []
    for i in range(df.shape[1]-1):
        ratio_list.append(py1_given_x(df.columns[i], "Label", df))

    return ratio_list
        


#===================================

def main():
    # load the training set
    df_train = load("netflix-train.csv", include_demographics=False)

    # compute model parameters (i.e. P(Y), P(X_i|Y))
    all_p_x_given_y = get_all_p_x_given_y("Label", df_train)

    x1y1 =[]
    for i in range(all_p_x_given_y.shape[0]):
        x1y1.append(round(all_p_x_given_y[i][1],8))
    
    #print(x1y1)

    p_y = get_p_y("Label", df_train)
    print(p_y)

    #calculate ratio for all movies
    ratio = three_movies(df_train)

    #find index of 3 movies with the highest ratio
    top3 = np.argsort(ratio)[-3:]
    for i in top3:
        print(df_train.columns[i])


    # load the test set
    df_test = load("netflix-test.csv", include_demographics=False)

    print(f"Training accuracy: {compute_accuracy(all_p_x_given_y, p_y, df_train)}")
    print(f"Test accuracy: {compute_accuracy(all_p_x_given_y, p_y, df_test)}")

if __name__ == "__main__":
    main()
