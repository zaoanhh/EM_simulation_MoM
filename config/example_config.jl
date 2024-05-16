if !isdefined(Main, :Drwatson)
    using DrWatson
end

freq = 5e9
grids_y = 1000
grids_x = 500
tag_size = 0.2
dis_target_rx = 3.0
doi_size_x = 5.0
doi_size_y = 10.0
tx_num = 1
rx_num = 3
dis_antenna = 0.03
dis_tx_rx = 3.0 + doi_size_x
N = 100
config = @strdict freq grids_x grids_y doi_size_x doi_size_y tx_num rx_num dis_antenna tag_size dis_tx_rx dis_target_rx N

