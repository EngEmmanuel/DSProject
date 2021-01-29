# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %%
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import time
from scipy.stats import pearsonr, spearmanr
from sklearn.model_selection import train_test_split, GridSearchCV, RandomizedSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.feature_selection import RFE
from sklearn.pipeline import Pipeline
from sklearn.metrics import r2_score
from sklearn.linear_model import LinearRegression, Ridge
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor

# Read in data

file_name = 'cah2_chembl_data_with_features-rdkit-Morgan_FP.csv' #'cah2_dude_data_plus_rdkit_descriptors.csv'
df = pd.read_csv(os.path.join('data',file_name))
df.drop_duplicates(subset=['molecule_chembl_id'], keep=False, inplace=True)
df.set_index('molecule_chembl_id', inplace=True)
df = df._get_numeric_data().dropna()


target = df.columns[0]
print(target)
y = df.pop(target).to_numpy()
X = df.to_numpy()
print(X.shape, y.shape, sep='\t')

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1, random_state=1)
X_train, X_val, y_train, y_val  = train_test_split(X_train, y_train, test_size=1/9, random_state=1)


# %%
def get_pipelines(features):
    p_lr=Pipeline([('sc1',StandardScaler()),
                         ('rfe1',RFE(estimator = DecisionTreeRegressor(), n_features_to_select=features[0],step= 0.56)),
                         ('lr_reg',LinearRegression())])

    p_gb=Pipeline([('sc2',StandardScaler()),
                         ('rfe2',RFE(estimator = DecisionTreeRegressor(), n_features_to_select=features[1],step= 0.56)),
                         ('gb_reg',GradientBoostingRegressor())])

    p_rfr=Pipeline([('sc3',StandardScaler()),
                         ('rfe3',RFE(estimator = DecisionTreeRegressor(), n_features_to_select=features[2],step= 0.56)),
                         ('rfr_reg',RandomForestRegressor())])

    return [p_lr, p_gb, p_rfr]


# %%

max_features = 200
featureRes = {}
for feature in range(max_features,max_features+1):

    features = [feature for i in range(3)]
    pipelines = get_pipelines(features)

    best_accuracy=0.0
    best_classifier=0
    best_pipeline=""

    # Dictionary of pipelines and classifier types for ease of reference
    pipe_dict = {0: 'Linear Regression', 1: 'Gradient Boost', 2: 'RandomForest'}

    # Fit the pipelines
    k = 0
    for pipe in pipelines:
    	pipe.fit(X_train, y_train)
    
    r2 = []
    for i,model in enumerate(pipelines):
        y_pred = model.predict(X_test)
        performance = model.score(X_test,y_test)
        r2.append(performance)
        print("{} Test Accuracy: {}".format(pipe_dict[i],performance))

    featureRes[str(feature)] = r2

    for i,model in enumerate(pipelines):
        if model.score(X_test,y_test)>best_accuracy:
            best_accuracy=model.score(X_test,y_test)
            best_pipeline=model
            best_classifier=i
    print('Classifier with best accuracy:{}'.format(pipe_dict[best_classifier]))

np.save('featureRes.npy',featureRes)

ylr_pred = pipelines[0].predict(X_test)
ygb_pred = pipelines[1].predict(X_test)
yrfr_pred = pipelines[2].predict(X_test)

np.savetxt("try_ygb_pred.csv", ygb_pred,delimiter= "\t")
np.savetxt("try_yrfr_pred.csv", yrfr_pred,delimiter= "\t")
np.savetxt("try_ylr_pred.csv", ylr_pred,delimiter= "\t")
np.savetxt("try_y_test.csv",y_test,delimiter= '\'t')

# %%
y = []
for k in range(len(pipelines)):
    temp = []
    for val in featureRes.values():
        temp.append(val[k])
    y.append(temp)
print(y)

for k in range(len(pipelines)):      
    plt.plot(featureRes.keys(), y[k],'r-')


#plt.close()


# %%

#optFeatures = [ y[i].index(max(y[i])) + 1 for i in range(3)]
#[p_lr, p_gb, p_rfr] = get_pipelines(optFeatures)

grid_param_rfr = {
                 "rfr_reg__n_estimators": [10, 100, 1000],
                 "rfr_reg__max_depth":[5,8,25,None],
                 "rfr_reg__min_samples_leaf":[1,2,5,15,100],
                 "rfr_reg__max_leaf_nodes": [2, 5,10]}

grid_param_gb =  {
                 "gb_reg__n_estimators": [10, 100, 1000],
                 "gb_reg__learning_rate":[0.5,0.8,0.25],
                 "gb_reg__min_samples_leaf":[1,2,5,15,100],
                 "gb_reg__max_leaf_nodes": [2,5,10]}

print('Starting randomised search')
gbsearch = RandomizedSearchCV(pipelines[1], grid_param_gb, n_iter= 50, cv = 3,verbose=2, n_jobs=-1)
rfrsearch = RandomizedSearchCV(pipelines[2], grid_param_rfr, n_iter= 50, cv = 3,verbose=2, n_jobs=-1)


best_gb_model = gbsearch.fit(X_train,y_train)
print('Finshed gradient boost fit')
best_rfr_model = rfrsearch.fit(X_train,y_train)

# Correlation plots


plt.scatter(y_test,ylr_pred)
plt.scatter(y_test,yrfr_pred)
plt.scatter(y_test,ygb_pred)

gb_pcorr, _ = pearsonr(y_test,ygb_pred)
gb_scorr, _ = spearmanr(y_test, ygb_pred)
rfr_pcorr, _ = spearmanr(y_test,yrfr_pred)
rfr_scorr, _ = pearsonr(y_test,yrfr_pred)
lr_pcorr, _ = spearmanr(y_test,ylr_pred)
lr_scorr, _ = pearsonr(y_test,ylr_pred)

print(gb_pcorr,gb_scorr,rfr_pcorr,rfr_scorr,lr_pcorr,lr_scorr)

print("The mean accuracy of the gb model is:",best_gb_model.score(X_test,y_test))
print("The mean accuracy of the rfr model is:",best_rfr_model.score(X_test,y_test))
print("The mean accuracy of the lr model is:",pipelines[0].score(X_test,y_test))
print(best_gb_model.best_estimator_,best_rfr_model.best_estimator_,sep = '\n\n')





# %%



