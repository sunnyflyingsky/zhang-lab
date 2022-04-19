#!/usr/bin/env python
# 
# Kui Xu, xukui.cs@gmail.com
# 2019-02-25
# ref smoothGrad

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.autograd import grad,Variable
import numpy as np
import os

class SmoothGrad(object):
    def __init__(self, model, device='cpu', only_seq=False, train=False, 
        x_stddev=0.015, y_stddev=0.015, nsamples=20, magnitude=2):
        self.model     = model
        self.device    = device
        self.train     = train
        self.only_seq  = only_seq
        self.x_stddev  = x_stddev
        self.y_stddev  = y_stddev
        self.nsamples  = nsamples
        self.magnitude = magnitude
        self.features  = model
        

    def get_gradients(self, X, y, pred_label=None):
        self.model.eval()
        self.model.zero_grad()
        X = X.to(self.device)
        y = y.to(self.device)
        #X.requires_grad=False

        y.requires_grad = True
        #output = self.model(z)
        #import pdb; pdb.set_trace()
        output = self.model(X, y)
        output = torch.sigmoid(output)
        output.backward()
        return y.grad

    def get_smooth_gradients(self, X, y, z=None):
        return self.__call__(X, y, z)
        
    def __call__(self, X, y, z=None):
        """[summary]
        
        Args:
            z ([type]): [description]
            y ([type]): [description]
            x_stddev (float, optional): [description]. Defaults to 0.15.
            t_stddev (float, optional): [description]. Defaults to 0.15.
            nsamples (int, optional):   [description]. Defaults to 20.
            magnitude (int, optional):  magnitude:0,1,2; 0: original gradient, 1: absolute value of the gradient,
                                        2: square value of the gradient. Defaults to 2.
        
        Returns:
            [type]: [description]
        """

        # 1. for sequece

        #x_stddev   = (self.x_stddev * (x.max()-x.min())).to(self.device).item()
        y_stddev = (self.y_stddev * (y.max() - y.min())).to(self.device).item()

        #total_x_grad = torch.zeros(x.shape).to(self.device)
        total_y_grad = torch.zeros(y.shape).to(self.device)
        #x_noise    = torch.zeros(x.shape).to(self.device)
        y_noise = torch.zeros(y.shape).to(self.device)
        """
        if not self.only_seq:
            # 2. for structure  
            #t = z[:,:,:,4:] #.data.cpu()
            t = z[:, :, :]  # .data.cpu()
            t_stddev = (self.t_stddev * (t.max()-t.min())).to(self.device).item()
            #t_total_grad = torch.zeros(t.shape)
            t_noise = torch.zeros(t.shape).to(self.device)
        """
        for i in range(self.nsamples):
           #x_plus_noise = x + x_noise.zero_().normal_(0, x_stddev)
            #y_plus_noise = y + y_noise.zero_().normal_(0, y_stddev)
            y_plus_noise = y
            """
            if self.only_seq:
                z_plus_noise = x_plus_noise
    
            else:
                t_plus_noise = t + t_noise.zero_().normal_(0, t_stddev)
                z_plus_noise = torch.cat((x_plus_noise, t_plus_noise), dim=3)
            """
            #print("z_plus_noise:",z_plus_noise.size())
            #grad = self.get_gradients(z_plus_noise, y)
            #import pdb; pdb.set_trace()
            y_grad = self.get_gradients(X, y_plus_noise, z)
            if self.magnitude == 1:
                #total_x_grad += torch.abs(x_grad)
                total_y_grad += torch.abs(y_grad)
            elif self.magnitude == 2:
                #total_x_grad += x_grad * x_grad
                total_y_grad += y_grad * y_grad

            # total_grad += grad * grad
        #total_x_grad /= self.nsamples
        total_y_grad /= self.nsamples
        return total_y_grad



    def get_batch_gradients(self, X, Y, Z=None):
        if Y is not None:
            #assert len(X) == len(Z), "The size of input {} and target {} are not matched.".format(len(X), len(Z))
            assert len(Y) == len(Z), "The size of input {} and target {} are not matched.".format(len(Y), len(Z))
        #d_g = torch.zeros_like(X)
        p_g = torch.zeros_like(Y)
        for i in range(Y.shape[0]):
            #x = X[i:i + 1]
            y = Y[i:i + 1]
            if Z is not None:
                z = Z[i:i+1]
            else:
                z = None
            p_g[i:i+1] = self.get_smooth_gradients(X,y,z)

        return p_g


def generate_saliency(model, x, y=None, smooth=False, nsamples=2, stddev=0.15, only_seq=False, train=False):
    saliency = SmoothGrad(model, only_seq, train)
    x_grad   = saliency.get_smooth_gradients(x, y, nsamples=nsamples, x_stddev=stddev, t_stddev=stddev)
    return x_grad



class GuidedBackpropReLU(torch.autograd.Function):

    def __init__(self, inplace=False):
        super(GuidedBackpropReLU, self).__init__()
        self.inplace = inplace

    def forward(self, input):
        pos_mask = (input > 0).type_as(input)
        output = torch.addcmul(
            torch.zeros(input.size()).type_as(input),
            input,
            pos_mask)
        self.save_for_backward(input, output)
        return output

    def backward(self, grad_output):
        input, output = self.saved_tensors

        pos_mask_1 = (input > 0).type_as(grad_output)
        pos_mask_2 = (grad_output > 0).type_as(grad_output)
        grad_input = torch.addcmul(
            torch.zeros(input.size()).type_as(input),
            torch.addcmul(
                torch.zeros(input.size()).type_as(input), grad_output, pos_mask_1),
                pos_mask_2)

        return grad_input

    def __repr__(self):
        inplace_str = ', inplace' if self.inplace else ''
        return self.__class__.__name__ + ' (' + inplace_str + ')'

class GuidedBackpropSmoothGrad(SmoothGrad):

    def __init__(self, model, device='cpu', only_seq=False, train=False, 
        x_stddev=0.15, t_stddev=0.15, nsamples=20, magnitude=2):
        super(GuidedBackpropSmoothGrad, self).__init__(
            model, device, only_seq, train, x_stddev, t_stddev, nsamples, magnitude)
        for idx, module in self.features._modules.items():
            if module.__class__.__name__ is 'ReLU':
                self.features._modules[idx] = GuidedBackpropReLU()

