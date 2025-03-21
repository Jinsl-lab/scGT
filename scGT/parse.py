from scGT import *

def parse_method(args, dataset, n, c, d, device):
    model=scGT(d, args.hidden_channels, c, num_layers=2, dropout=0.0,
                num_heads=1, use_bn=True, nb_random_features=35,
                use_gumbel=True, use_residual=True, use_act=False, use_jk=False,
                nb_gumbel_sample=15, rb_order=2, rb_trans='sigmoid').to(device)
    return model

