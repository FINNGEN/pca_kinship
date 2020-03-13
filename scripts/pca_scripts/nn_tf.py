from utils import basic_iterator,np,mapcount,return_header,tmp_bash,merge_files,make_sure_path_exists,return_header,identify_separator,write_fam_samplelist
import pandas as pd
import os,shlex,subprocess

import torch
import torch.nn as nn
import torch.nn.functional as F
device = torch.device("cpu")

def load_tensor_data(args):

    tag = 'FIN'
    eigenvec_file = args.pca_outlier_path+ tag + ".eigenvec"
    tmp_data = pd.read_csv(eigenvec_file,identify_separator(eigenvec_file))
    fin_inputs = torch.Tensor(tmp_data.values)
    fin_targets = torch.ones(fin_inputs.shape[0],dtype = torch.float32)
    
    tag = 'NON_FIN'
    eigenvec_file = args.pca_outlier_path+ tag + ".eigenvec"
    tmp_data = pd.read_csv(eigenvec_file,identify_separator(eigenvec_file))
    non_fin_inputs = torch.Tensor(tmp_data.values)
    non_fin_targets = torch.zeros(non_fin_inputs.shape[0],dtype = torch.float32)

    train_inputs= torch.cat((fin_inputs, non_fin_inputs), 0)
    train_targets = torch.cat((fin_targets, non_fin_targets), 0)
    
    
    tag = 'FINNGEN'
    eigenvec_file = args.pca_outlier_path+ tag + ".eigenvec"
    tmp_data = pd.read_csv(eigenvec_file,identify_separator(eigenvec_file))
    finngen_inputs = torch.Tensor(tmp_data.values)
    finngen_targets = torch.ones(finngen_inputs.shape[0],dtype = torch.float32)

 
    # Normalize inputs to zero mean and unit variance
    mean = train_inputs.mean(dim=0)
    std = train_inputs.std(dim=0)
    scaler = lambda x: (x - mean.to(x.device)) / std.to(x.device)
    
    train_inputs = scaler(train_inputs)
    finngen_inputs = scaler(finngen_inputs)

    print(train_inputs.shape)
    print(train_inputs)
    print(train_targets.shape)
    print(train_targets)
    
    return train_inputs,train_targets,finngen_inputs,finngen_targets

class MLP(nn.Module):
    def __init__(self, sizes, activation_fn=F.tanh):
        """Multilayer perceptron with an arbitrary number of layers.
        
        Args:
          sizes (list): Number of units in each layer including the input and the output layer:
                         [n_inputs, n_units_in_hidden_layer1, ..., n_units_in_hidden_layerN, n_outputs]
          activation_fn (callable): An element-wise function used in every layer except in the last one.
          out_fn (callable): An element-wise function used in the last layer. By default a lambda func that returns the output of Linear

        """
        super(MLP, self).__init__()
        self.hidden = nn.ModuleList()
        *hidden_sizes,out_size = sizes
        for k in range(len(hidden_sizes)-1):
            self.hidden.append(nn.Linear(hidden_sizes[k], hidden_sizes[k+1]))
        # Output layer
        self.out = nn.Linear(hidden_sizes[-1], out_size)
        self.activation = activation_fn
        
    def forward(self, x):
        # YOUR CODE HERE
        for layer in self.hidden:
            x = self.activation(layer(x))
        return torch.sigmoid(self.out(x))

def compute_accuracy(model, inputs, targets):
    with torch.no_grad():
        inputs, targets = inputs.to(device), targets.to(device)
        outputs = (model.forward(inputs) > 0.5).float()
        accuracy = (outputs == targets).sum().float() / targets.numel()
        return accuracy

    
def model_run(args):
    
    # Create the model
    model = MLP([3, 100,50, 10, 1])
    model.to(device)
    print(model)

    train_inputs,train_targets,finngen_inputs,finngen_targets = load_tensor_data(args)

    # set gradients to zero
    iterations = 2000
    optimizer = torch.optim.Adam(model.parameters(),lr = 0.01)
    print(f"Accuracy at step 0: {compute_accuracy(model,finngen_inputs,finngen_targets)}")
    for i in range(1,1+iterations):
        # calculate bce and backdrop   
        model.zero_grad()
        out = model(train_inputs)
        loss = F.binary_cross_entropy(out,train_targets)
        loss.backward()
        optimizer.step()
        if not i % 200 or i == iterations:
            print(f"Accuracy at step {i}: {compute_accuracy(model,train_inputs,train_targets)}")  
    print(f"Accuracy of test targets: {compute_accuracy(model,finngen_inputs,finngen_targets)}")
    
def merged_pca(args):
    #######
    # PCA #
    #######
    args.pca_output_file = os.path.join(args.pca_outlier_path,args.name)
    if not os.path.isfile( args.pca_output_file+ '.eigenval') or args.force:
        args.force = True 
        #individuals that need to be removed
        cmd = f'plink2 --bfile {args.merged_plink_file} --read-freq  {args.merged_plink_file}.afreq  --pca {args.pca_components} approx biallelic-var-wts --threads {args.cpus}  -out {pca_output_file}'
        print(cmd)
        subprocess.call(shlex.split(cmd))




def build_train_data(args):
    '''
    Writes pop info for outlier detection
    '''
    tg_pop = os.path.join(args.data_path,'20130606_sample_info.txt')

    columns = range(5,5+args.pc_filter)
    columns = list(map(str,columns))
    cut_columns =  ','.join(columns)  

    tag = 'FIN'
    eigenvec_file = args.pca_outlier_path+ tag + ".eigenvec"
    if not os.path.isfile(eigenvec_file):
        # PROJECT FINNS
        fin_fam = os.path.join(args.pca_outlier_path,'FIN.fam')
        cmd = f"cat {tg_pop} | grep -w FIN | cut -f 1,2 > {fin_fam}"
        tmp_bash(cmd)
    
        cmd = f'plink2 --bfile {args.new_tg} --keep {fin_fam} --score {args.pca_output_file+".eigenvec.var"} 2 3 header-read no-mean-imputation  --score-col-nums {columns[0]}-{columns[-1]} --out {args.pca_outlier_path+tag}'
        print(cmd)
        subprocess.call(shlex.split(cmd))
        cmd = f'cat {args.pca_outlier_path+tag + ".sscore"}  | cut -f {cut_columns} >  {eigenvec_file}'
        tmp_bash(cmd)

    tag = 'NON_FIN'
    eigenvec_file = args.pca_outlier_path+ tag + ".eigenvec"
    if not os.path.isfile(eigenvec_file):
        # PROJECT FINNS
        fin_fam = os.path.join(args.pca_outlier_path,'FIN.fam')
        cmd = f"cat {tg_pop} | grep -wv FIN | cut -f 1,2 > {fin_fam}"
        tmp_bash(cmd)
    
        cmd = f'plink2 --bfile {args.new_tg} --remove {fin_fam} --score {args.pca_output_file+".eigenvec.var"} 2 3 header-read no-mean-imputation  --score-col-nums {columns[0]}-{columns[-1]} --out {args.pca_outlier_path+tag}'
        print(cmd)
        subprocess.call(shlex.split(cmd))
        cmd = f'cat {args.pca_outlier_path+tag + ".sscore"}  | cut -f{cut_columns} >  {eigenvec_file}'
        tmp_bash(cmd)

    tag = 'FINNGEN'
    eigenvec_file = args.pca_outlier_path+ tag + ".eigenvec"
    if not os.path.isfile(eigenvec_file):
        cmd = f'plink2 --bfile {args.merged_plink_file} --keep {args.sample_fam} --score {args.pca_output_file+".eigenvec.var"} 2 3 header-read no-mean-imputation  --score-col-nums {columns[0]}-{columns[-1]} --out {args.pca_outlier_path+tag}'
        print(cmd)
        subprocess.call(shlex.split(cmd))
        cmd = f'cat {args.pca_outlier_path+tag + ".sscore"}  | cut -f {cut_columns} >  {args.pca_outlier_path+ tag + ".eigenvec"}'
        tmp_bash(cmd)
    
