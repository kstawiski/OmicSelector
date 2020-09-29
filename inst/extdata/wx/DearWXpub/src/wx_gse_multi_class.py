import os
import argparse
import numpy as np
import pandas as pd
import _pickle as cPickle #for pyton3.x
from sklearn.utils import shuffle
from sklearn import cross_validation
from keras.utils import to_categorical
from wx_hyperparam import WxHyperParameter
from wx_core import wx_slp, wx_mlp, classifier_LOOCV
from tqdm import tqdm

def StandByRow(data_frame, unused_cols=[]):
    unused_cols_list = []
    for col_ in unused_cols:
        unused_cols_list.append(data_frame[col_])
    data_frame = data_frame.drop(unused_cols,axis=1)

    data_frame = data_frame.astype(float)
    data_frame.fillna(0, inplace=True)
    data_frame = data_frame.apply(lambda x: ((x-x.mean())/x.std()), axis=1)

    for n,col_ in enumerate(unused_cols):
        data_frame[col_] = unused_cols_list[n]

    return data_frame


def make_data_frame_GSE105127(raw_file = './', sel_class=[], norm_flag = True):
    df = pd.read_csv(raw_file,sep='\t')
    df = df.drop(['chr','start','end','strand','Length','Copies','Annotation/Divergence'],axis=1)
    cols = df.columns

    cols_cv = []
    cols_iz = []
    cols_pp = []
    cnt_cv=0
    cnt_iz=0
    cnt_pp=0
    for n,name in enumerate(cols):
        clf = name.split('_')
        if len(clf) > 3:
            clf = clf[2]
            if clf == sel_class[0]:
                cnt_cv+=1
                cols_cv.append(name)
            if clf == sel_class[1]:
                cnt_iz+=1            
                cols_iz.append(name)
            if clf == sel_class[2]:
                cnt_pp+=1            
                cols_pp.append(name)                
    print(cnt_cv, cnt_iz, cnt_pp)

    fix_col_name = ['Transcript_ID']
    df = df.set_index(fix_col_name[0]).T#set cell column to index
    df[fix_col_name[0]] = df.index
    df = df.reset_index(drop=True)

    if norm_flag:
        df = StandByRow(df, unused_cols=fix_col_name)

    f_names = df.columns.tolist()
    f_names.remove('Transcript_ID')

    df_cv = df[df['Transcript_ID'].isin(cols_cv)]
    df_iz = df[df['Transcript_ID'].isin(cols_iz)]
    df_pp = df[df['Transcript_ID'].isin(cols_pp)]

    # df_con = df_con.drop('Gene', axis=1)
    # df_ad = df_ad.drop('Gene', axis=1)

    #add label
    cv_label = np.zeros(df_cv.shape[0])
    iz_label = np.ones(df_iz.shape[0])
    pp_label = np.empty(df_pp.shape[0])
    pp_label.fill(2)

    df_cv['label'] = cv_label
    df_iz['label'] = iz_label
    df_pp['label'] = pp_label        

    df = pd.concat([df_cv, df_iz, df_pp],axis=0)
    df = df.rename(columns={'Transcript_ID':'id'})


    return df, f_names


def make_data_frame_GSE112057(raw_folder = './', sel_class=[], norm_flag = True):
    data_file = 'GSE112057_Normalized_dataset.txt'
    anno_file = 'GSE112057_clinic.xlsx'

    df_data = pd.read_csv(raw_folder+data_file,sep='\t')
    f_names = df_data.gene.values

    df_anno = pd.read_excel(raw_folder+anno_file)

    fix_col_name = ['gene']
    df_data = df_data.set_index(fix_col_name[0]).T#set cell column to index
    df_data[fix_col_name[0]] = df_data.index
    df_data = df_data.reset_index(drop=True)
    df_data_ids = df_data.gene.values
    df_data_labels = np.empty(len(df_data_ids))
    df_data_labels.fill(-1)


    if norm_flag:
        df_data = StandByRow(df_data, unused_cols=fix_col_name)


    for n,class_name in enumerate(sel_class):
        class_ids = df_anno[df_anno['disease state (diagnosis)'] == class_name].ID.values
        df_data_labels[np.where(np.isin(df_data_ids, class_ids))] = n
        
    df_data['label'] = df_data_labels
    df_data = df_data[df_data.label != -1]
    df_data = df_data.rename(columns={"gene": "id"})


    return df_data, f_names


def load_norm_feature_set(df, validation_ratio, RANDOM_STATE, num_cls):
    def get_val_and_train_at_label(label):
        df_label = df[df.label == label]
        df_label = df_label.drop(['label'],axis=1)

        label_values = df_label.values.tolist()
        label_values = shuffle(label_values, random_state = RANDOM_STATE)    

        label_th = int(float(len(label_values)) * validation_ratio)

        val_x = label_values[:label_th]
        train_x = label_values[label_th:]       


        train_y = []
        val_y = []
        for i in range(0,label_th):
            val_y.append(label)
        for i in range(label_th,len(label_values)):
            train_y.append(label)        

        return val_x, val_y, train_x, train_y 

    all_val_x = []
    all_val_y = []
    all_train_x = []
    all_train_y = []    
    for label_at in range(0, num_cls):
        v_x,v_y,t_x,t_y = get_val_and_train_at_label(label_at)
        all_val_x = all_val_x + v_x
        all_val_y = all_val_y + v_y
        all_train_x = all_train_x + t_x
        all_train_y = all_train_y + t_y           

    #shuffle
    all_train_x, all_train_y = shuffle(all_train_x, all_train_y, random_state = RANDOM_STATE)
    all_val_x, all_val_y = shuffle(all_val_x, all_val_y, random_state = RANDOM_STATE)

    all_train_y = to_categorical(all_train_y, num_cls)
    all_val_y = to_categorical(all_val_y, num_cls)

    return np.asarray(all_train_x), np.asarray(all_train_y), np.asarray(all_val_x), np.asarray(all_val_y) 
    

def wx_feature_selection(df='', gene_names='', n_sel = 14, val_ratio = 0.2, iter=1000, epochs=30, learning_ratio=0.001, batch_size=32, 
                        verbose=False, model_type='MLP', num_cls=2):

    feature_num = len(gene_names)
    all_weight = np.zeros(feature_num)    
    all_count = np.ones(feature_num)
    for i in range(0, iter):
        train_x, train_y, val_x, val_y = load_norm_feature_set(df, val_ratio, i, num_cls)
        print(i, 'train : ',train_x.shape, 'val : ',val_x.shape)
        hp = WxHyperParameter(epochs=epochs, learning_ratio=learning_ratio, batch_size=batch_size, verbose=verbose)
        if model_type == 'MLP':
            sel_idx, sel_weight, val_acc = wx_mlp(train_x, train_y, val_x, val_y, n_selection=min(n_sel*100, feature_num), hyper_param=hp, num_cls=num_cls)
        if model_type == 'SLP':
            sel_idx, sel_weight, val_acc = wx_slp(train_x, train_y, val_x, val_y, n_selection=min(n_sel*100, feature_num), hyper_param=hp, num_cls=num_cls)
        for j in range(0,min(n_sel*100, feature_num)):
            all_weight[sel_idx[j]] += sel_weight[j]
            all_count[sel_idx[j]] += 1
    all_weight = all_weight / all_count
    sort_index = np.argsort(all_weight)[::-1]    
    sel_index = sort_index[:n_sel]
    sel_weight =  all_weight[sel_index]
    gene_names = np.asarray(gene_names)
    sel_genes = gene_names[sel_index]

    return sel_index, sel_genes, sel_weight

def evaluation_LOOCV(df, sel_genes, gene_names, method_clf='xgb', verbose=False, norm_method='no',num_cls=2):
    RANDOM_STATE = 1
    VAL_RATIO = 0.2
    # VAL_RATIO = 0.01
    
    def do_LOOCV(all_x, all_y):
        loo = cross_validation.LeaveOneOut(len(all_x))
        tot_cnt = np.zeros(num_cls)
        hit_cnt = np.zeros(num_cls)

        cancer_prob = []
        labels = []
        cnt=0
        for train_index, test_index in tqdm(loo):
            train_val_x, test_x = all_x[train_index], all_x[test_index]
            train_val_y, test_y = all_y[train_index], all_y[test_index]

            train_val_x, train_val_y = shuffle(train_val_x, train_val_y, random_state=RANDOM_STATE)
            n_trn = len(train_val_x)
            n_dev = int(n_trn*VAL_RATIO)
            n_trn = n_trn - n_dev
            train_x = train_val_x[0:n_trn]
            train_y = train_val_y[0:n_trn]
            val_x = train_val_x[n_trn:]
            val_y = train_val_y[n_trn:]     

            prob = classifier_LOOCV(train_x, train_y, val_x, val_y, test_x, test_y, method_clf=method_clf, verbose=verbose, num_cls=num_cls)

            for i in range(0,num_cls):
                if test_y[0] == i:
                    tot_cnt[i]+=1
                    if np.argmax(prob) == i:
                        hit_cnt[i]+=1

        return (tot_cnt,hit_cnt)

    type_acc = 0

    sel_genes.append('label')
    if norm_method == 'no':
        df = df[sel_genes]       
    if norm_method == 'after':
        df = df[sel_genes]               
        df = StandByRow(df,['label'])
    if norm_method == 'before':        
        df = StandByRow(df,['label'])
        df = df[sel_genes]               


    data_label = df['label'].values.astype(int) #tumor or normal label
    #print(data_label.shape, sum(data_label))

    df = df.drop(['label'],axis=1)
    data_x = df.values #seleted features

    return do_LOOCV(data_x, data_label)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--gse', type=int, default=105127, help='gse number')
    parser.add_argument('--gpu', type=int, default=0, help='gpu num')    
    parser.add_argument('--n_sel', type=int, default=30, help='sel gene num')
    parser.add_argument('--eval_gene', type=int, default=30, help='eval gene num')    
    args, unparsed = parser.parse_known_args()

    gse_number = args.gse
    gpu_number = args.gpu
    n_sel = args.n_sel
    eval_gene_number = args.eval_gene

    os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_number)

    if gse_number == 112057:
        sel_class = ['Crohn Disease', 'Oligoarticular JIA', 'Polyarticular JIA', 'Systemic JIA', 'Ulcerative Colitis', 'Control']
        df, f_names = make_data_frame_GSE112057('./GSE_DATA/', sel_class = sel_class, norm_flag=False)
    if gse_number == 105127:
        sel_class = ['CV','IZ','PP']
        df, f_names = make_data_frame_GSE105127('./GSE_DATA/Total_hg38.txt', sel_class = sel_class, norm_flag=True)

    def get_before_df(label):
        df_label = df[df.label == label]
        th = df_label.shape[0]
        return df_label[:int(th/2)]
    def get_after_df(label):
        df_label = df[df.label == label]
        th = df_label.shape[0]
        return df_label[int(th/2):th]

    concat_list = []
    for i in range(0,len(sel_class)):
        concat_list.append(get_before_df(i))
    df_feature_select = pd.concat(concat_list,axis=0)

    concat_list = []
    for i in range(0,len(sel_class)):
        concat_list.append(get_after_df(i))    
    df_eval = pd.concat(concat_list,axis=0)
    print(df_feature_select.shape)
    print(df_eval.shape)

    df_feature_select = df_feature_select.drop('id', axis=1)
    df_eval = df_eval.drop('id', axis=1)

    if True:
        sel_idx, sel_genes, sel_weight = wx_feature_selection(df=df_feature_select, gene_names=f_names, n_sel = n_sel, val_ratio = 0.2, iter=10, epochs=30,
                                    learning_ratio=0.001, batch_size=16, verbose=False, model_type='SLP', num_cls=len(sel_class))
        print ('\nSingle Layer WX')
        print ('selected feature index:',sel_idx)
        print ('selected feature genes:',sel_genes)
        print ('selected feature weights:',sel_weight)

        cPickle.dump(sel_genes,open('gse'+str(gse_number)+'_sel_genes.cpickle','wb'))

    sel_genes = cPickle.load(open('gse'+str(gse_number)+'_sel_genes.cpickle','rb'))
    sel_genes = sel_genes[:eval_gene_number]

    tot_cnt, hit_cnt = evaluation_LOOCV(df_eval, sel_genes.tolist(), f_names, method_clf='xgb', verbose=False, norm_method='no', num_cls=len(sel_class))

    for i in range(0,len(tot_cnt)):
        print(sel_class[i], int(tot_cnt[i]), int(hit_cnt[i]), float(hit_cnt[i]/tot_cnt[i]*100))
    print(float(np.sum(hit_cnt)/np.sum(tot_cnt)*100))