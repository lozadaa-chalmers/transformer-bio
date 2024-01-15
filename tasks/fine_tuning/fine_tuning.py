import loralib as lora
from torch import nn

from torch.utils.data import DataLoader
from typing import Literal


def fine_tune(model: nn.Module,
              training_loader: DataLoader = None,
              validation_loader: DataLoader = None,
              optimizer:
              strategy: Literal['lora'] = None,
              lora_r: int = 8,
              d_model: int = None,
              bias: Literal['none', 'all', 'lora_only'] = 'none'):
    used_model = model
    if strategy == 'lora':
        lora.mark_only_lora_as_trainable(model=used_model,
                                         bias=bias)
        for data, label in training_loader:


        qkv_proj = lora.MergedLinear(in_features=d_model,
                                     out_features=3 * d_model,
                                     r=lora_r,
                                     enable_lora=[True, False, True])
