#ifndef ENGINE_NET
#define ENGINE_NET

#include "array.hpp"
#include "global.hpp"

//--------------------------------------
struct nnet
{
    array1d<float> input_mean;
    array1d<float> input_std;
    array1d<float> output_mean;
    array1d<float> output_std;
    std::vector<array2d<float>> weights;
    std::vector<array1d<float>> biases;
};

void nnet_load(nnet& nn, const char* filename)
{
    FILE* f = fopen(filename, "rb");
    assert(f != NULL);
    
    array1d_read(nn.input_mean, f);
    array1d_read(nn.input_std, f);
    array1d_read(nn.output_mean, f);
    array1d_read(nn.output_std, f);
    
    int count;
    fread(&count, sizeof(int), 1, f);
    
    nn.weights.resize(count);
    nn.biases.resize(count);
    
    for(int i = 0; i < count; i++)
    {
        array2d_read(nn.weights[i], f);
        array1d_read(nn.biases[i], f);
    }
    
    fclose(f);
}

//--------------------------------------

static inline void nnet_layer_normalize(
    slice1d<float> output,
    const slice1d<float> mean,
    const slice1d<float> std)
{
    for(int i = 0; i < output.size; i++)
    {
        output(i) = (output(i) - mean(i)) / std(i);
    }
}

static inline void nnet_layer_denormalize(
    slice1d<float> output,
    const slice1d<float> mean,
    const slice1d<float> std)
{
    for(int i = 0; i < output.size; i++)
    {
        output(i) = output(i) * std(i) + mean(i);
    }
}

static inline void nnet_layer_linear(
    slice1d<float> output,
    const slice1d<float> input,
    const slice2d<float> weights,
    const slice1d<float> biases)
{
    for(int j = 0; j < output.size; j++)
    {
        output(j) = biases(j);
    }
    
    for(int i = 0; i < input.size; i++)
    {
        if (input(i) != 0.0f)
        {
            for(int j = 0; j < output.size; j++)
            {
                output(j) += input(i) * weights(i, j);
            }
        }
    }
}

static inline void nnet_layer_relu(slice1d<float> output)
{
    for(int i = 0; i < output.size; i++)
    {
        output(i) = std::max(output(i), 0.f);
    }
}

static inline void nnet_layer_elu(slice1d<float> output)
{
    for(int i = 0; i < output.size; i++)
    {
        if(output(i) < 0.f)
        {
            output(i) = exp(output(i)) - 1.f;
        }
    }
}

//--------------------------------------
struct nnet_evaluation
{
    std::vector<array1d<float>> layers;
    
    // Resize for a given network
    void resize(const nnet& nn)
    {
        layers.resize(nn.weights.size() + 1);
        layers.front().resize(nn.weights.front().rows);
      
        for(int i = 0; i < nn.weights.size(); i++)
        {
            layers[i + 1].resize(nn.weights[i].cols);            
        }
    }
};

void nnet_evaluate(
    nnet_evaluation& evaluation,
    const nnet& nn)
{
    nnet_layer_normalize(
        evaluation.layers.front(),
        nn.input_mean,
        nn.input_std);    
  
    for(int i = 0; i < nn.weights.size(); i++)
    {
        nnet_layer_linear(
            evaluation.layers[i + 1],
            evaluation.layers[i],
            nn.weights[i],
            nn.biases[i]);
        
        // No relu for final layer
        if (i != nn.weights.size() - 1)
        {
            nnet_layer_elu(evaluation.layers[i + 1]);
        }
    }
    
    nnet_layer_denormalize(
        evaluation.layers.back(),
        nn.output_mean,
        nn.output_std);
}

void compressor_evaluate(
    slice1d<float> rotations,
    slice1d<float> positions,
    float t,
    bool ifcontroll,
    nnet_evaluation& evaluation,
    const nnet& nn)
{
    if(!ifcontroll)
    {
        rotations(0) = cos(deg_to_rad(2.2 * t)); 
        rotations(1) = circulate_float(0.1 * t, -180.f, 180.f) / 180.f;
        rotations(2) = 0.f;
    }

    slice1d<float> input_layer = evaluation.layers.front();
    slice1d<float> output_layer = evaluation.layers.back();

    for(int i = 0; i < rotations.size; i++)
    {
        input_layer(i) = rotations(i);
    }
    
    nnet_evaluate(evaluation, nn);
    
    for(int i = 0; i < positions.size; i++)
    {
        positions(i) = output_layer(i);
    }
}

void decompressor_evaluate(
    slice1d<float> positions,
    slice1d<float> rotations,
    nnet_evaluation& evaluation,
    const nnet& nn)
{
    slice1d<float> input_layer = evaluation.layers.front();
    slice1d<float> output_layer = evaluation.layers.back();
  
    for(int i = 0; i < positions.size; i++)
    {
        input_layer(i) = positions(i);
    }
    
    nnet_evaluate(evaluation, nn);
    
    for(int i = 0; i < rotations.size; i++)
    {
        rotations(i) = output_layer(i);
    }
}

#endif