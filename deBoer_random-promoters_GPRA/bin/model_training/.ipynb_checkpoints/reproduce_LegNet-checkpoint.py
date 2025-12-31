import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader

class SEBlock(nn.Module):
    def __init__(self, channel, reduction=16):
        super(SEBlock, self).__init__()
        self.avg_pool = nn.AdaptiveAvgPool1d(1)
        self.fc = nn.Sequential(
            nn.Linear(channel, channel // reduction, bias=False),
            nn.ReLU(inplace=True),
            nn.Linear(channel // reduction, channel, bias=False),
            nn.Sigmoid()
        )

    def forward(self, x):
        b, c, _ = x.size()
        y = self.avg_pool(x).view(b, c)
        y = self.fc(y).view(b, c, 1)
        return x * y.expand_as(x)

class EfficientNetBlock(nn.Module):
    def __init__(self, in_channels, out_channels, groups):
        super(EfficientNetBlock, self).__init__()
        self.conv = nn.Conv1d(in_channels, out_channels, kernel_size=3, stride=1, padding=1, groups=groups)
        self.bn = nn.BatchNorm1d(out_channels)
        self.silu = nn.SiLU()
        self.se = SEBlock(out_channels)
        self.residual = nn.Conv1d(in_channels, out_channels, kernel_size=1) if in_channels != out_channels else nn.Identity()

    def forward(self, x):
        residual = self.residual(x)
        x = self.conv(x)
        x = self.bn(x)
        x = self.silu(x)
        x = self.se(x)
        return torch.cat([x, residual], dim=1)  # Channel-wise concatenation

class LegNet(nn.Module):
    def __init__(self):
        super(LegNet, self).__init__()
        self.stem = nn.Sequential(
            nn.Conv1d(6, 256, kernel_size=7, stride=1, padding=3, bias=False),
            nn.BatchNorm1d(256),
            nn.SiLU()
        )
        channels = [256, 128, 128, 64, 64, 64, 64]
        self.blocks = nn.Sequential(*[EfficientNetBlock(channels[i], channels[i+1], groups=1) for i in range(6)])
        self.final_conv = nn.Conv1d(channels[-1], 10, kernel_size=1)  # Adjust the number of output channels as needed
        self.pool = nn.AdaptiveAvgPool1d(1)
        self.softmax = nn.Softmax(dim=1)

    def forward(self, x):
        x = self.stem(x)
        x = self.blocks(x)
        x = self.final_conv(x)
        x = self.pool(x)
        x = torch.flatten(x, 1)
        x = self.softmax(x)
        return x