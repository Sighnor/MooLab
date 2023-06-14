import sys
import struct

import numpy as np

import matplotlib.pyplot as plt

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.tensorboard import SummaryWriter

from train_common import load_features, load_frames, save_network

# Networks

class Decompressor(nn.Module):

    def __init__(self, input_size, output_size, hidden_size=64):
        super(Decompressor, self).__init__()
        
        self.linear0 = nn.Linear(input_size, hidden_size)
        self.linear1 = nn.Linear(hidden_size, hidden_size)
        self.linear2 = nn.Linear(hidden_size, hidden_size)
        self.linear3 = nn.Linear(hidden_size, hidden_size)
        self.linear4 = nn.Linear(hidden_size, output_size)

    def forward(self, x):
        x = F.elu(self.linear0(x))
        x = F.elu(self.linear1(x))
        x = F.elu(self.linear2(x))
        x = F.elu(self.linear3(x))
        return self.linear4(x)

# Training procedure

if __name__ == '__main__':
    
    # Load data
    
    Feature = load_features('resources/features.bin')['features'].copy().astype(np.float32)
    Frame = load_frames('resources/frames.bin')['frames'].copy().astype(np.float32)
    
    size = Feature.shape[0]
    nfeatures = Feature.shape[1]
    nframes = Frame.shape[1]
    
    # Parameters
    
    seed = 1234
    batchsize = 4
    lr = 0.001
    niter = 500000
    
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.set_num_threads(1)
    
    # Compute means/stds
    
    Feature_scale = Feature.std()
    Feature_noise_std = Feature.std(axis=0)
    
    decompressor_mean_out = torch.as_tensor(np.hstack([
        Frame.mean(axis=0).ravel(),
    ]).astype(np.float32))
    
    decompressor_std_out = torch.as_tensor(np.hstack([
        Frame.std(axis=0).ravel(),
    ]).astype(np.float32))
    
    decompressor_mean_in = torch.as_tensor(np.hstack([
        Feature.mean(axis=0).ravel(),
    ]).astype(np.float32))
    
    decompressor_std_in = torch.as_tensor(np.hstack([
        Feature_scale.repeat(nfeatures),
    ]).astype(np.float32))
    
    # Make networks
    
    network_decompressor = Decompressor(nfeatures, nframes)
    
    # Function to generate test predictions

    def generate_predictions():
        
        with torch.no_grad():
            
            # Get slice of database for first clip
            
            start = 0
            stop = size
            
            nsigma = np.random.uniform(size=[stop-start, 1]).astype(np.float32)
            noise = np.random.normal(size=[stop-start, nfeatures]).astype(np.float32)
            Featurehat = Feature[start:stop] + 0.05 * Feature_noise_std * noise
            
            # Find nearest
            
            Framegnd = torch.as_tensor(Frame[start:stop])
            Featurehat = torch.as_tensor(Featurehat)
            
            # Decompress
            
            output = (network_decompressor((Featurehat - decompressor_mean_in) / decompressor_std_in) *
                decompressor_std_out + decompressor_mean_out)
            
            Frametil = output[:,:]
            
            # Write frames
            
            lmin, lmax = Framegnd.cpu().numpy().min(), Framegnd.cpu().numpy().max()
            
            fig, axs = plt.subplots(nframes, sharex=True, figsize=(12, 2*nframes))
            for i in range(nframes):
                axs.plot(Framegnd[:2500:5,i].cpu().numpy(), marker='.', linestyle='None')
                axs.plot(Frametil[:2500:5,i].cpu().numpy(), marker='.', linestyle='None')
                axs.set_ylim(lmin, lmax)
            plt.tight_layout()
            
            try:
                plt.savefig('resources/feature_frame.png')
            except IOError as e:
                print(e)

            plt.close()

    # Train
    
    writer = SummaryWriter()

    optimizer = torch.optim.AdamW(
        network_decompressor.parameters(), 
        lr=lr,
        amsgrad=True,
        weight_decay=0.001)
        
    scheduler = torch.optim.lr_scheduler.ExponentialLR(optimizer, gamma=0.99)
    
    rolling_loss = None
    
    sys.stdout.write('\n')
    
    for i in range(niter):
    
        optimizer.zero_grad()
        
        # Extract batch
        
        samples = np.random.randint(0, size, size=[batchsize])
        nsigma = np.random.uniform(size=[batchsize, 1]).astype(np.float32)
        noise = np.random.normal(size=[batchsize, nfeatures]).astype(np.float32)
        
        Featurehat = Feature[samples] + 0.05 * Feature_noise_std * noise
        
        # Find frame
        
        Framegnd = Frame[samples]

        Featurehat = torch.as_tensor(Featurehat)
        Framegnd = torch.as_tensor(Framegnd)
        
        # Decompressor
        
        output = (network_decompressor((Featurehat - decompressor_mean_in) / decompressor_std_in) *
            decompressor_std_out + decompressor_mean_out)
        
        Frametil = output[:,:]
        
        # Compute Losses
        
        loss_frameval = torch.mean(torch.abs(Framegnd - Frametil))
        loss = loss_frameval
        
        # Backprop
        
        loss.backward()

        optimizer.step()
    
        # Logging
        
        writer.add_scalar('resources/loss', loss.item(), i)
        
        writer.add_scalars('resources/loss_terms', {
            'frameval': loss_frameval.item(),
        }, i)
        
        if rolling_loss is None:
            rolling_loss = loss.item()
        else:
            rolling_loss = rolling_loss * 0.99 + loss.item() * 0.01
        
        if i % 10 == 0:
            sys.stdout.write('\rIter: %7i Loss: %5.3f' % (i, rolling_loss))
        
        if i % 1000 == 0:
            generate_predictions()
            save_network('resources/feature_frame.bin', [
                network_decompressor.linear0, 
                network_decompressor.linear1, 
                network_decompressor.linear2, 
                network_decompressor.linear3,
                network_decompressor.linear4],
                decompressor_mean_in,
                decompressor_std_in,
                decompressor_mean_out,
                decompressor_std_out)
            
        if i % 1000 == 0:
            scheduler.step()
            