import tensorflow as tf
from keras.models import Model
from keras.layers import Input, Dense
from keras import backend as K
from keras import optimizers,applications, callbacks
from keras.callbacks import ModelCheckpoint
from keras.callbacks import LearningRateScheduler
import numpy as np
from wx_hyperparam import WxHyperParameter
import xgboost as xgb
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
import functools
import time

#set default global hyper paramerters
wx_hyperparam = WxHyperParameter(learning_ratio=0.001)

def timeit(func):
    @functools.wraps(func)
    def newfunc(*args, **kwargs):
        startTime = time.time()
        ret = func(*args, **kwargs)
        elapsedTime = time.time() - startTime
        # print('function [{}] finished in {} ms'.format(
        #     func.__name__, int(elapsedTime * 1000)))
        print('\nfunction [{}] finished in {} s'.format(
            func.__name__, float(elapsedTime)))
        return ret
    return newfunc

def cw_ann_model(x_train, y_train, x_val, y_val, hyper_param=wx_hyperparam, hidden_layer_size=128, num_cls=2):
    input_dim = len(x_train[0])
    inputs = Input((input_dim,))
    hidden = Dense(hidden_layer_size)(inputs)
    fc_out = Dense(num_cls,  activation='softmax')(hidden)
    model = Model(input=inputs, output=fc_out)
    # model.summary()

    #build a optimizer
    sgd = optimizers.SGD(lr=hyper_param.learning_ratio, decay=hyper_param.weight_decay, momentum=hyper_param.momentum, nesterov=True)
    model.compile(loss='categorical_crossentropy', optimizer=sgd, metrics=['accuracy'])        

    #call backs
    def step_decay(epoch):
        exp_num = int(epoch/10)+1       
        return float(hyper_param.learning_ratio/(10 ** exp_num))

    best_model_path="../slp_cw_ann_weights_best"+".hdf5"
    save_best_model = ModelCheckpoint(best_model_path, monitor="val_loss", verbose=hyper_param.verbose, save_best_only=True, mode='min')
    #save_best_model = ModelCheckpoint(best_model_path, monitor="val_acc", verbose=1, save_best_only=True, mode='max')
    change_lr = LearningRateScheduler(step_decay)                                

    #run train
    history = model.fit(x_train, y_train, validation_data=(x_val,y_val), 
                epochs=hyper_param.epochs, batch_size=hyper_param.batch_size, shuffle=True, callbacks=[save_best_model, change_lr], verbose=hyper_param.verbose)

    #load best model
    model.load_weights(best_model_path)

    return model

@timeit
def connection_weight(x_train, y_train, x_val, y_val, n_selection=100, hidden_layer_size=128, hyper_param=wx_hyperparam, num_cls=2):
    input_dim = len(x_train[0])

    # make model and do train
    model = cw_ann_model(x_train, y_train, x_val, y_val, hyper_param=hyper_param, hidden_layer_size=hidden_layer_size, num_cls=num_cls)

    #load weights
    weights = model.get_weights()

    #get feature importance using connection weight algo (Olden 2004)
    wt_ih = weights[0]#.transpose() #input-hidden weights
    wt_ho = weights[1]#.transpose() #hidden-out weights
    dot_wt = wt_ih * wt_ho
    sum_wt = np.sum(dot_wt,axis=1)

    selected_idx = np.argsort(sum_wt)[::-1][0:n_selection]
    selected_weights = sum_wt[selected_idx]

    #get evaluation acc from best model
    loss, val_acc = model.evaluate(x_val, y_val)

    K.clear_session()

    return selected_idx, selected_weights, val_acc

def naive_SLP_model(x_train, y_train, x_val, y_val, hyper_param=wx_hyperparam, num_cls=2):
    input_dim = len(x_train[0])
    inputs = Input((input_dim,))
    #fc_out = Dense(2,  kernel_initializer='zeros', bias_initializer='zeros', activation='softmax')(inputs)
    fc_out = Dense(num_cls,  activation='softmax')(inputs)
    model = Model(input=inputs, output=fc_out)
    # model.summary()

    #build a optimizer
    sgd = optimizers.SGD(lr=hyper_param.learning_ratio, decay=hyper_param.weight_decay, momentum=hyper_param.momentum, nesterov=True)
    model.compile(loss='categorical_crossentropy', optimizer=sgd, metrics=['accuracy'])        

    #call backs
    def step_decay(epoch):
        exp_num = int(epoch/10)+1       
        return float(hyper_param.learning_ratio/(10 ** exp_num))

    best_model_path="../slp_wx_weights_best"+".hdf5"
    save_best_model = ModelCheckpoint(best_model_path, monitor="val_loss", verbose=hyper_param.verbose, save_best_only=True, mode='min')
    #save_best_model = ModelCheckpoint(best_model_path, monitor="val_acc", verbose=1, save_best_only=True, mode='max')
    change_lr = LearningRateScheduler(step_decay)                                

    #run train
    history = model.fit(x_train, y_train, validation_data=(x_val,y_val), 
                epochs=hyper_param.epochs, batch_size=hyper_param.batch_size, shuffle=True, callbacks=[save_best_model, change_lr], verbose=hyper_param.verbose)

    #load best model
    model.load_weights(best_model_path)

    return model

def classifier_LOOCV(x_train, y_train, x_val, y_val, x_test, y_test, method_clf='xgb', verbose=False, num_cls=2):
    if method_clf=='xgb':
        if num_cls == 2:
            clf = xgb.XGBClassifier(seed=1, objective='binary:logistic')
            clf.fit(x_train, y_train, eval_set=[(x_val, y_val)], verbose=verbose, eval_metric='logloss', early_stopping_rounds=100)
            pred_prob = clf.predict_proba(x_test)
            return pred_prob[0][1]
        else:
            clf = xgb.XGBClassifier(seed=1, objective='multi:softprob')
            clf.fit(x_train, y_train, eval_set=[(x_val, y_val)], verbose=verbose, eval_metric='mlogloss', early_stopping_rounds=100)
            pred_prob = clf.predict_proba(x_test)
            return pred_prob[0]
        

    if method_clf=='svm':
        # clf = SVC(kernel = 'linear')
        if num_cls == 2:
            clf = SVC(kernel='rbf', probability=True, C=1.0, degree=3, verbose=verbose, random_state=0)
            #print(x_train.shape, y_train.shape, x_test.shape)
            clf.fit(x_train,y_train)
            pred_prob = clf.predict_proba(x_test)
            return pred_prob[0][1]
        else:
            clf = SVC(kernel='rbf', probability=True, C=1.0, degree=3, verbose=verbose, random_state=0)
            #print(x_train.shape, y_train.shape, x_test.shape)
            clf.fit(x_train,y_train)
            pred_prob = clf.predict_proba(x_test)            
            return pred_prob[0]
    

@timeit
def wx_slp(x_train, y_train, x_val, y_val, n_selection=100, hyper_param=wx_hyperparam, num_cls=2):
    if num_cls < 2:
        return

    # sess = tf.Session()
    # K.set_session(sess)

    input_dim = len(x_train[0])

    # make model and do train
    model = naive_SLP_model(x_train, y_train, x_val, y_val, hyper_param=hyper_param, num_cls=num_cls)

    #load weights
    weights = model.get_weights()

    #cacul WX scores
    num_data = {}
    running_avg={}
    tot_avg={}
    Wt = weights[0].transpose() #all weights of model
    Wb = weights[1].transpose() #all bias of model
    for i in range(num_cls):
        tot_avg[i] = np.zeros(input_dim) # avg of input data for each output class
        num_data[i] = 0.
    for i in range(len(x_train)):
        c = y_train[i].argmax()
        x = x_train[i]
        tot_avg[c] = tot_avg[c] + x
        num_data[c] = num_data[c] + 1
    for i in range(num_cls):
        tot_avg[i] = tot_avg[i] / num_data[i]

    #for general multi class problems
    wx_mul = []
    for i in range(0,num_cls):
        wx_mul_at_class = []
        for j in range(0,num_cls):
            wx_mul_at_class.append(tot_avg[i] * Wt[j])
        wx_mul.append(wx_mul_at_class)
    wx_mul = np.asarray(wx_mul)

    wx_abs = np.zeros(Wt.shape[1])
    for n in range(0, Wt.shape[1]):
        for i in range(0,num_cls):
            for j in range(0,num_cls):
                if i != j:
                    wx_abs[n] += np.abs(wx_mul[i][i][n] - wx_mul[i][j][n])

    selected_idx = np.argsort(wx_abs)[::-1][0:n_selection]
    selected_weights = wx_abs[selected_idx]

    #get evaluation acc from best model
    loss, val_acc = model.evaluate(x_val, y_val)

    K.clear_session()

    return selected_idx, selected_weights, val_acc


def sum_fan_in(xi, input_x, layer_num, index, wt, output_class_idx):
    #wx = ux*uw
    # print('call ', index)
    if index == layer_num - 1:#backprop output layer
        cur_x = sum_fan_in(xi, input_x, layer_num, index-1, wt, output_class_idx)
        cur_w = wt[index][output_class_idx]
        cur_wx = cur_x * cur_w
        ret = np.sum(cur_wx)

    elif index == 0:#handle input layer
        cur_x = input_x[xi]
        cur_w = wt[index]
        cur_wx = []
        for i in range(0,len(wt[index])):#loop for hidden units
            cur_wx.append(cur_x * wt[index][i][xi]) 
        ret = np.asarray(cur_wx)

    else:#normal hiddenlayer backprop
        cur_x = sum_fan_in(xi, input_x, layer_num, index-1, wt, output_class_idx)
        cur_wx = []
        for i in range(0,len(wt[index])):#loop for hidden units 
            local_wx = cur_x * wt[index][i]
            local_sum = np.sum(local_wx)
            cur_wx.append(local_sum)
        ret = np.asarray(cur_wx)            

    return ret

def cal_class_wx_mlp(input_avg, wt, wb, input_class_idx, output_class_idx):
    layer_num = len(wt)
    num_feature = len(input_avg[0])
    wx = []
    for i in range(0, num_feature):
        wx.append( sum_fan_in(i, input_avg[input_class_idx], layer_num, layer_num - 1, wt, output_class_idx) )

    #print('mlp wx_done ', input_class_idx, output_class_idx)
    return np.asarray(wx)

@timeit
def wx_mlp(x_train, y_train, x_val, y_val, n_selection=100, hyper_param=wx_hyperparam, num_cls=2):
    if num_cls < 2:
        return

    #sess = tf.Session()
    #K.set_session(sess)    
    
    #build a NN model
    input_dim = len(x_train[0])
    num_hidden_layer = hyper_param.num_hidden_layer
    num_h_unit = hyper_param.num_h_unit
    inputs = Input((input_dim,))
    hidden_1 = Dense(units=num_h_unit)(inputs)
    hidden_2 = Dense(units=int(num_h_unit/2))(hidden_1)    
    fc_out = Dense(num_cls,  kernel_initializer='zeros', bias_initializer='zeros', activation='softmax')(hidden_2)
    model = Model(input=inputs, output=fc_out)

    #build a optimizer
    sgd = optimizers.SGD(lr=hyper_param.learning_ratio, decay=hyper_param.weight_decay, momentum=hyper_param.momentum, nesterov=True)
    model.compile(loss='categorical_crossentropy', optimizer=sgd, metrics=['accuracy'])    
    # model.summary()

    #call backs
    def step_decay(epoch):
        exp_num = int(epoch/10)+1       
        return float(hyper_param.learning_ratio/(10 ** exp_num))

    best_model_path="./mlp_wx_weights_best"+".hdf5"
    save_best_model = ModelCheckpoint(best_model_path, monitor="val_loss", verbose=hyper_param.verbose, save_best_only=True, mode='min')
    #save_best_model = ModelCheckpoint(best_model_path, monitor="val_acc", verbose=1, save_best_only=True, mode='max')
    change_lr = LearningRateScheduler(step_decay)                                

    #run train
    history = model.fit(x_train, y_train, validation_data=(x_val,y_val), 
                epochs=hyper_param.epochs, batch_size=hyper_param.batch_size, shuffle=True, callbacks=[save_best_model, change_lr], verbose=hyper_param.verbose)

    #load best model
    model.load_weights(best_model_path)

    #load weights
    weights = model.get_weights()

    #cacul WX scores
    num_data = {}
    running_avg={}
    tot_avg={}
    # load weights and bias
    wt = {}
    wb = {}
    for i in range(0,num_hidden_layer+1):
        wt[i] = weights[i*2].transpose() 
        wb[i] = weights[i*2+1].transpose()
    # make avg for input data
    for i in range(num_cls):
        tot_avg[i] = np.zeros(input_dim) # avg of input data for each output class
        num_data[i] = 0.
    for i in range(len(x_train)):
        c = y_train[i].argmax()
        x = x_train[i]
        tot_avg[c] = tot_avg[c] + x
        num_data[c] = num_data[c] + 1
    for i in range(num_cls):
        tot_avg[i] = tot_avg[i] / num_data[i]        

    #for general multi class problems
    wx_mul = []
    for i in range(0,num_cls):
        wx_mul_at_class = []
        for j in range(0,num_cls):
            #wx_mul_at_class.append(tot_avg[i] * Wt[j])
            print('Cal mlp wx : input class, weight class = ',i,j)
            wx_mul_at_class.append( cal_class_wx_mlp(tot_avg, wt, wb, i, j) )
        wx_mul.append(wx_mul_at_class)
    wx_mul = np.asarray(wx_mul)

    wx_abs = np.zeros(input_dim)
    for n in range(0, input_dim):
        for i in range(0,num_cls):
            for j in range(0,num_cls):
                if i != j:
                    wx_abs[n] += np.abs(wx_mul[i][i][n] - wx_mul[i][j][n])


    selected_idx = np.argsort(wx_abs)[::-1][0:n_selection]
    selected_weights = wx_abs[selected_idx]

    #get evaluation acc from best model
    loss, val_acc = model.evaluate(x_val, y_val)

    K.clear_session()

    return selected_idx, selected_weights, val_acc    
