import os
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from keras.utils import to_categorical
from wx_hyperparam import WxHyperParameter
from wx_core import wx_slp, wx_mlp, connection_weight


def get_sample_data(num_cls=2):
    train_num = 100
    test_num = 100
    input_dim = 20000
    num_cls = num_cls
    if num_cls < 2:
        return

    x_train = np.random.random((train_num, input_dim))
    y_train = to_categorical(np.random.randint(num_cls, size=(train_num, 1)), num_classes=num_cls)

    x_test = np.random.random((test_num, input_dim))
    y_test = to_categorical(np.random.randint(num_cls, size=(test_num, 1)), num_classes=num_cls)

    return x_train, y_train, x_test, y_test


def StandByCol(data_frame, unused_cols=[]):
    unused_cols_list = []
    for col_ in unused_cols:
        unused_cols_list.append(data_frame[col_])
    data_frame = data_frame.drop(unused_cols,axis=1)

    data_frame = data_frame.astype(float)
    data_frame.fillna(0, inplace=True)
    #data_frame = data_frame.apply(lambda x: ((x-x.mean())/x.std()), axis=0)

    for n,col_ in enumerate(unused_cols):
        data_frame[col_] = unused_cols_list[n]

    return data_frame
    

if __name__ == '__main__':
    os.environ["CUDA_VISIBLE_DEVICES"] = "0"    

    #random values
    if False:
        num_cls = 5
        x_train, y_train, x_val, y_val = get_sample_data(num_cls = num_cls)

    #read from actual datas
    df_data = pd.read_csv('test-data.csv')
    f_names = df_data['fnames'].values
    df_data = df_data.drop(['fnames'],axis=1)
    id_names = df_data.columns.values
    
    #z-scoring
    df_data = StandByCol(df_data,unused_cols=[])
    df_data = df_data.T
    x_val = df_data.values    

    df_anno = pd.read_csv('test-anno.csv')
    anno_ids = df_anno.id.values
    anno_class = df_anno.label.values
    class_type = df_anno.label.unique()

    n_cls = len(class_type)
    y_val = []
    for id_ in id_names:
        idx = np.where(anno_ids == id_)
        y_val.append(np.where(class_type == anno_class[idx])[0][0])
    y_val = np.asarray(y_val)
    y_val = to_categorical(y_val, num_classes=n_cls)

    print('TEST samples names : ', id_names)    
    print('TEST classes : ' , class_type)    

    #read from actual datas
    df_data = pd.read_csv('train-data.csv')
    f_names = df_data['fnames'].values
    df_data = df_data.drop(['fnames'],axis=1)
    id_names = df_data.columns.values
    
    #z-scoring
    df_data = StandByCol(df_data,unused_cols=[])
    df_data = df_data.T
    x_train = df_data.values    

    df_anno = pd.read_csv('train-anno.csv')
    anno_ids = df_anno.id.values
    anno_class = df_anno.label.values
    class_type = df_anno.label.unique()

    n_cls = len(class_type)
    y_train = []
    for id_ in id_names:
        idx = np.where(anno_ids == id_)
        y_train.append(np.where(class_type == anno_class[idx])[0][0])
    y_train = np.asarray(y_train)
    y_train = to_categorical(y_train, num_classes=n_cls)

    print('TRAIN samples names : ', id_names)    
    print('TRAIN classes : ' , class_type)    

    #split to train and val
    #x_train, x_val, y_train, y_val = train_test_split(x_all, y_all, test_size=0.2, random_state=1)

    hp = WxHyperParameter(epochs=3080, learning_ratio=0.01, batch_size=100, num_hidden_layer = 2, num_h_unit = 16, verbose=True)
    sel_idx, sel_weight, val_acc = wx_slp(x_train, y_train, x_val, y_val, n_selection=10, hyper_param=hp, num_cls=n_cls)

    print ('\nSingle Layer WX')
    print ('selected feature names:',f_names[sel_idx])
    print ('selected feature index:',sel_idx)
    print ('selected feature weight:',sel_weight)
    print ('evaluation accuracy:',val_acc)
    print ('\n\n\n\n')
    
    file1 = open("wyniki.txt","w") 
    file1.write(np.array2string(f_names[sel_idx],separator=','))
    file1.close()
    
    np.save("wyniki.npy", f_names[sel_idx])
