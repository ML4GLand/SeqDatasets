class EffBlock(nn.Module):
    def __init__(self, 
        in_ch, 
        ks, 
        resize_factor,
        filter_per_group,
        activation, 
        out_ch=None,
        se_reduction=None,
        se_type="simple",
        inner_dim_calculation="in"
    ):
        super().__init__()
        self.in_ch = in_ch
        self.out_ch = self.in_ch if out_ch is None else out_ch
        self.resize_factor = resize_factor
        self.se_reduction = resize_factor if se_reduction is None else se_reduction
        self.ks = ks
        self.inner_dim_calculation = inner_dim_calculation

        '''
        `in` refers to the original method of EfficientNetV2 to set the dimensionality of the EfficientNetV2-like block
        `out` is the mode used in the original LegNet approach

        This parameter slighly changes the mechanism of channel number calculation 
        which can be seen in the figure above (C, channel number is highlighted in red).
        '''
        if inner_dim_calculation == "out":
            self.inner_dim = self.out_ch * self.resize_factor
        elif inner_dim_calculation == "in":
            self.inner_dim = self.in_ch * self.resize_factor
        else:
            raise Exception(f"Wrong inner_dim_calculation: {inner_dim_calculation}")
            
        self.filter_per_group = filter_per_group

        se_constructor = SELayer

        block = nn.Sequential(
            nn.Conv1d(
                in_channels=self.in_ch,
                out_channels=self.inner_dim,
                kernel_size=1,
                padding='same',
                bias=False
           ),
           nn.BatchNorm1d(self.inner_dim),
           activation(),
           nn.Conv1d(
                in_channels=self.inner_dim,
                out_channels=self.inner_dim,
                kernel_size=ks,
                groups=self.inner_dim // self.filter_per_group,
                padding='same',
                bias=False
           ),
           nn.BatchNorm1d(self.inner_dim),
           activation(),
           se_constructor(
               self.in_ch, 
              self.inner_dim,
              reduction=self.se_reduction
           ),
           nn.Conv1d(
                in_channels=self.inner_dim,
                out_channels=self.in_ch,
                kernel_size=1,
                padding='same',
                bias=False
           ),
           nn.BatchNorm1d(self.in_ch),
           activation(),
        )
        self.block = block
    
    def forward(self, x):
        return self.block(x)
    
    
class MappingBlock(nn.Module):
    def __init__(self, in_ch, out_ch, activation):
        super().__init__()
        self.block = nn.Sequential(
            nn.Conv1d(
                in_channels=in_ch,
                out_channels=out_ch,
                kernel_size=1,
                padding='same',
           ),
           activation()
        )
        
    def forward(self, x):
        return self.block(x)
    
    
class ResidualConcat(nn.Module):
    def __init__(self, fn):
        super().__init__()
        self.fn = fn

    def forward(self, x, **kwargs):
        return torch.concat([self.fn(x, **kwargs), x], dim=1)