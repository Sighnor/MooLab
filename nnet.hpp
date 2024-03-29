#ifndef MOOLAB_NET
#define MOOLAB_NET

#include "array.hpp"
#include "global.hpp"

// Taken from https://github.com/orangeduck/Motion-Matching
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

void nnet_evaluate(
    slice1d<float> input,
    slice1d<float> output,
    nnet_evaluation& evaluation,
    const nnet& nn)
{
    slice1d<float> input_layer = evaluation.layers.front();
    slice1d<float> output_layer = evaluation.layers.back();

    assert(input.size == input_layer.size && output.size == output_layer.size);

    for(int i = 0; i < input.size; i++)
    {
        input_layer(i) = input(i);
    }
    
    nnet_evaluate(evaluation, nn);
    
    for(int i = 0; i < output.size; i++)
    {
        output(i) = output_layer(i);
    }
}

void nnet_evaluate(
    slice1d<float> input,
    slice1d<quat> output,
    nnet_evaluation& evaluation,
    const nnet& nn)
{
    slice1d<float> input_layer = evaluation.layers.front();
    slice1d<float> output_layer = evaluation.layers.back();

    assert(input.size == input_layer.size && (4 * output.size) == output_layer.size);

    for(int i = 0; i < input.size; i++)
    {
        input_layer(i) = input(i);
    }
    
    nnet_evaluate(evaluation, nn);
    
    for(int i = 0; i < output.size; i++)
    {
        output(i).w = output_layer(4 * i);
        output(i).x = output_layer(4 * i + 1);
        output(i).y = output_layer(4 * i + 2);
        output(i).z = output_layer(4 * i + 3);
    }
}

void nnet_evaluate(
    slice1d<quat> input,
    slice1d<float> output,
    nnet_evaluation& evaluation,
    const nnet& nn)
{
    slice1d<float> input_layer = evaluation.layers.front();
    slice1d<float> output_layer = evaluation.layers.back();

    assert((4 * input.size) == input_layer.size && output.size == output_layer.size);

    for(int i = 0; i < input.size; i++)
    {
        input_layer(4 * i) = input(i).w;
        input_layer(4 * i + 1) = input(i).x;
        input_layer(4 * i + 2) = input(i).y;
        input_layer(4 * i + 3) = input(i).z;
    }
    
    nnet_evaluate(evaluation, nn);
    
    for(int i = 0; i < output.size; i++)
    {
        output(i) = output_layer(i);
    }
}

#endif
