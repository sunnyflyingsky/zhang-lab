import os

import torch.nn.functional as F
import torch

def set(drug_encoding=None,
        target_encoding=None,
        result_folder="smartnet/smartnet_vscode/result/",
        input_dim_drug=2586,
        input_dim_protein=6,
        hidden_dim_drug=256,
        hidden_dim_protein=256,
        cls_hidden_dims=[1024, 1024, 512],
        mlp_hidden_dims_drug=[1024, 256, 64],
        mlp_hidden_dims_target=[1024, 256, 64],
        batch_size=256,
        train_epoch=10,
        #test_every_X_epoch=20,
        LR=1e-4,
        decay=0,
        transformer_emb_size_drug=128,
        transformer_intermediate_size_drug=512,
        transformer_num_attention_heads_drug=8,
        transformer_n_layer_drug=8,
        transformer_emb_size_target=64,
        transformer_intermediate_size_target=256,
        transformer_num_attention_heads_target=4,
        transformer_n_layer_target=2,
        transformer_dropout_rate=0.1,
        transformer_attention_probs_dropout=0.1,
        transformer_hidden_dropout_rate=0.1,
        mpnn_hidden_size=50,
        mpnn_depth=3,
        cnn_drug_filters=[32, 64, 96],
        cnn_drug_kernels=[4, 6, 8],
        cnn_target_filters=[32, 64, 96],
        cnn_target_kernels=[4, 8, 12],
        rnn_Use_GRU_LSTM_drug='GRU',
        rnn_drug_hid_dim=64,
        rnn_drug_n_layers=2,
        rnn_drug_bidirectional=True,
        rnn_Use_GRU_LSTM_target='GRU',
        rnn_target_hid_dim=64,
        rnn_target_n_layers=2,
        rnn_target_bidirectional=True,
        num_workers=0,
        inter=0,

#=============================================================================
        gnn_hid_dim_drug = 64,
        gnn_num_layers = 3,
        gnn_activation = F.relu,
        in_feats = 74,
#=============================================================================
        neuralfp_max_degree=10,
        neuralfp_predictor_hid_dim=128,
        neuralfp_predictor_activation=torch.tanh,
        attentivefp_num_timesteps=2,
#=============================================================================
        gat_num_heads=10,
        gat_feat_drops=0.2,
        gat_attn_drops=0.2,
#=============================================================================


        ):
    base_config = {'input_dim_drug': input_dim_drug,
                   'input_dim_protein': input_dim_protein,
                   'hidden_dim_drug': hidden_dim_drug,  # hidden dim of drug
                   'hidden_dim_protein': hidden_dim_protein,  # hidden dim of protein
                   'cls_hidden_dims': cls_hidden_dims,  # decoder classifier dim 1
                   'batch_size': batch_size,
                   'train_epoch': train_epoch,
                   # 'test_every_X_epoch': test_every_X_epoch,
                   'LR': LR,
                   'drug_encoding': drug_encoding,
                   'target_encoding': target_encoding,
                   'result_folder': result_folder,
                   'num_workers': num_workers,
                   'inter': inter,
                   }
    if not os.path.exists(base_config['result_folder']):
        os.makedirs(base_config['result_folder'])
    """
    if drug_encoding == 'Pubchem':
        base_config['input_dim_drug'] = 881
        base_config['mlp_hidden_dims_drug'] = mlp_hidden_dims_drug  # MLP classifier dim 1

    elif drug_encoding == 'Daylight':
        base_config['input_dim_drug'] = 2048
        base_config['mlp_hidden_dims_drug'] = mlp_hidden_dims_drug  # MLP classifier dim 1
    """

    if drug_encoding == 'CNN':
        base_config['cnn_drug_filters'] = cnn_drug_filters
        base_config['cnn_drug_kernels'] = cnn_drug_kernels



    elif drug_encoding == 'RNN':
        base_config['rnn_Use_GRU_LSTM_drug'] = rnn_Use_GRU_LSTM_drug
        base_config['rnn_drug_hid_dim'] = rnn_drug_hid_dim
        base_config['rnn_drug_n_layers'] = rnn_drug_n_layers
        base_config['rnn_drug_bidirectional'] = rnn_drug_bidirectional
        base_config['cnn_drug_filters'] = cnn_drug_filters
        base_config['cnn_drug_kernels'] = cnn_drug_kernels

    elif drug_encoding == 'Transformer':
        #base_config['input_dim_drug'] = 2586
        base_config['input_dim_drug'] = 64
        base_config['transformer_emb_size_drug'] = transformer_emb_size_drug
        base_config['transformer_num_attention_heads_drug'] = transformer_num_attention_heads_drug
        base_config['transformer_intermediate_size_drug'] = transformer_intermediate_size_drug
        base_config['transformer_n_layer_drug'] = transformer_n_layer_drug
        base_config['transformer_dropout_rate'] = transformer_dropout_rate
        base_config['transformer_attention_probs_dropout'] = transformer_attention_probs_dropout
        base_config['transformer_hidden_dropout_rate'] = transformer_hidden_dropout_rate
        base_config['hidden_dim_drug'] = transformer_emb_size_drug


    elif drug_encoding == 'MPNN':
        base_config['hidden_dim_drug'] = hidden_dim_drug
        base_config['batch_size'] = batch_size
        base_config['mpnn_hidden_size'] = mpnn_hidden_size
        base_config['mpnn_depth'] = mpnn_depth





# =============================================================================
    elif drug_encoding == 'DGL_GCN':
        base_config['gnn_hid_dim_drug'] = gnn_hid_dim_drug
        base_config['gnn_num_layers'] = gnn_num_layers
        base_config['gnn_activation'] = gnn_activation
        base_config['in_feats'] = in_feats
        # =============================================================================
    elif drug_encoding == 'DGL_GAT':
        base_config['gnn_hid_dim_drug'] = gnn_hid_dim_drug
        base_config['gnn_num_layers'] = gnn_num_layers
        base_config['gnn_activation'] = gnn_activation
        base_config['gat_num_heads'] = gat_num_heads
        base_config['gat_feat_drops'] = gat_feat_drops
        base_config['gat_attn_drops'] = gat_attn_drops
        base_config['in_feats'] = in_feats

    # =============================================================================
    elif drug_encoding == 'DGL_NeuralFP':
        base_config['gnn_hid_dim_drug'] = gnn_hid_dim_drug
        base_config['neuralfp_max_degree'] = neuralfp_max_degree
        base_config['gnn_activation'] = gnn_activation
        base_config['neuralfp_predictor_hid_dim'] = neuralfp_predictor_hid_dim
        base_config['gnn_num_layers'] = gnn_num_layers
        base_config['neuralfp_predictor_activation'] = neuralfp_predictor_activation
        base_config['in_feats'] = in_feats
    # =============================================================================
    elif drug_encoding == 'DGL_AttentiveFP':
        base_config['gnn_hid_dim_drug'] = gnn_hid_dim_drug
        base_config['gnn_num_layers'] = gnn_num_layers
        base_config['attentivefp_num_timesteps'] = attentivefp_num_timesteps

    elif drug_encoding == 'DGL_GIN_AttrMasking':
        ## loaded pretrained model specifications
        pass
    elif drug_encoding == 'DGL_GIN_ContextPred':
        ## loaded pretrained model specifications
        pass


    # =============================================================================

    # raise NotImplementedError
    elif drug_encoding is None:
        pass
    else:
        raise AttributeError("Please use the correct drug encoding available!")

    if target_encoding == 'AAC':
        base_config['mlp_hidden_dims_target'] = mlp_hidden_dims_target  # MLP classifier dim 1
    elif target_encoding == 'PseudoAAC':
        base_config['input_dim_protein'] = 30
        base_config['mlp_hidden_dims_target'] = mlp_hidden_dims_target  # MLP classifier dim 1
    elif target_encoding == 'Conjoint_triad':
        base_config['input_dim_protein'] = 343
        base_config['mlp_hidden_dims_target'] = mlp_hidden_dims_target  # MLP classifier dim 1
    elif target_encoding == 'Quasi-seq':
        base_config['input_dim_protein'] = 100
        base_config['mlp_hidden_dims_target'] = mlp_hidden_dims_target  # MLP classifier dim 1

    elif target_encoding == 'CNN':
        base_config['cnn_target_filters'] = cnn_target_filters
        base_config['cnn_target_kernels'] = cnn_target_kernels
    elif target_encoding == 'Prism':
        base_config['cnn_target_filters'] = cnn_target_filters
        base_config['cnn_target_kernels'] = cnn_target_kernels

    elif target_encoding == 'RNN':
        base_config['rnn_Use_GRU_LSTM_target'] = rnn_Use_GRU_LSTM_target
        base_config['rnn_target_hid_dim'] = rnn_target_hid_dim
        base_config['rnn_target_n_layers'] = rnn_target_n_layers
        base_config['rnn_target_bidirectional'] = rnn_target_bidirectional
        base_config['cnn_target_filters'] = cnn_target_filters
        base_config['cnn_target_kernels'] = cnn_target_kernels
    elif target_encoding == 'Transformer':
        #base_config['input_dim_protein'] = 4114
        base_config['input_dim_protein'] = 6
        base_config['transformer_emb_size_target'] = transformer_emb_size_target
        base_config['transformer_num_attention_heads_target'] = transformer_num_attention_heads_target
        base_config['transformer_intermediate_size_target'] = transformer_intermediate_size_target
        base_config['transformer_n_layer_target'] = transformer_n_layer_target
        base_config['transformer_dropout_rate'] = transformer_dropout_rate
        base_config['transformer_attention_probs_dropout'] = transformer_attention_probs_dropout
        base_config['transformer_hidden_dropout_rate'] = transformer_hidden_dropout_rate
        base_config['hidden_dim_protein'] = transformer_emb_size_target
    elif target_encoding is None:
        pass
    else:
        raise AttributeError("Please use the correct protein encoding available!")

    return base_config
