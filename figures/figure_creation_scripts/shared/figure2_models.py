"""Feedforward model used by the supplied DeepFlux checkpoints."""

from torch import nn


class MLPLayer(nn.Module):
    def __init__(self, in_dim, out_dim, activation, use_batch_norm, dropout_rate):
        super().__init__()
        self.linear = nn.Linear(in_dim, out_dim)
        self.bn = nn.BatchNorm1d(out_dim) if use_batch_norm else None
        self.activation = activation
        self.dropout = nn.Dropout(dropout_rate) if dropout_rate > 0 else nn.Identity()

    def forward(self, x):
        x = self.linear(x)
        if self.bn is not None:
            x = self.bn(x)
        return self.dropout(self.activation(x))


class DeepFlux(nn.Module):
    def __init__(self, parameters):
        super().__init__()
        self.use_residual = parameters["use_residual"]
        activation = {"relu": nn.ReLU, "tanh": nn.Tanh, "sigmoid": nn.Sigmoid, "leaky_relu": nn.LeakyReLU}[parameters["activation_function"]]()
        layers = []
        incoming = parameters["input_dim"]
        outgoing = parameters["dim_size"]
        for _ in range(parameters["num_layers"]):
            layers.append(MLPLayer(incoming, outgoing, activation, parameters["use_batch_norm"], parameters["dropout_rate"]))
            incoming = outgoing
            outgoing = max(1, int(outgoing * parameters["dim_shrinking_scale"]))
        self.shared_layers = nn.ModuleList(layers)
        self.head = nn.Linear(incoming, parameters["output_dim"])

    def forward(self, x):
        for layer in self.shared_layers:
            new = layer(x)
            x = x + new if self.use_residual and x.shape == new.shape else new
        return self.head(x)
