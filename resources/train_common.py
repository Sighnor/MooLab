import struct
import numpy as np
import torch        

def load_frames(filename):

    with open(filename, 'rb') as f:
        
        size, nframes = struct.unpack('II', f.read(8))
        size = 2500
        frames = np.frombuffer(f.read(size * nframes * 4), dtype = np.float32, count = size * nframes).reshape([size, nframes])
        
    return {
        'frames': frames,
    }

def load_features(filename):

    with open(filename, 'rb') as f:
        
        size, nfeatures = struct.unpack('II', f.read(8))
        size = 2500
        features = np.frombuffer(f.read(size * nfeatures * 4), dtype = np.float32, count = size * nfeatures).reshape([size, nfeatures])
        
    return {
        'features': features,
    }

    
def save_network(filename, layers, mean_in, std_in, mean_out, std_out):
    
    with torch.no_grad():
        
        with open(filename, 'wb') as f:
            f.write(struct.pack('I', mean_in.shape[0]) + mean_in.cpu().numpy().astype(np.float32).ravel().tobytes())
            f.write(struct.pack('I', std_in.shape[0]) + std_in.cpu().numpy().astype(np.float32).ravel().tobytes())
            f.write(struct.pack('I', mean_out.shape[0]) + mean_out.cpu().numpy().astype(np.float32).ravel().tobytes())
            f.write(struct.pack('I', std_out.shape[0]) + std_out.cpu().numpy().astype(np.float32).ravel().tobytes())
            f.write(struct.pack('I', len(layers)))
            for layer in layers:
                f.write(struct.pack('II', *layer.weight.T.shape) + layer.weight.T.cpu().numpy().astype(np.float32).ravel().tobytes())
                f.write(struct.pack('I', *layer.bias.shape) + layer.bias.cpu().numpy().astype(np.float32).ravel().tobytes())
