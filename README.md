# EM_simulation_MoM
Simulate the propagation of two-dimensional TM waves using the julia language and the method of moments.


## 文件说明
- `src/`目录下是常用函数
- `notebook/`目录下是jupyter notebook文件
- `config/`目录下是仿真参数配置文件
- `data/`目录下是仿真结果数据

## 依赖安装
```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

## 使用说明
示例代码在`notebook/example.ipynb`中，可以通过jupyter notebook打开查看。
1. 首先导入相关函数
    ```julia
    using DrWatson
    @quickactivate "EIT"
    using Revise, ProgressMeter, Distributions, Random
    includet(srcdir("funcs.jl"))
    import ..Simu
    using ..MoMForward: ConstantParameter
    import ..MoMForward as mf
    ```

    在这里我们使用了 [DrWatson](https://github.com/JuliaDynamics/DrWatson.jl) 这个包来做项目中代码的管理。它有几个便于使用的点：1.通过 `@quickactivate` 来激活项目环境；2.通过 `srcdir` 来指定源代码目录等函数来对对应的文件夹做了映射，这样无论待执行脚本处于什么路径中，都可以通过这种方式来引入对应的文件。3. `produce_or_load` 函数可以用来判断是否需要重新生成数据，如果数据已经存在则直接加载，否则生成新的数据。详情查看它的文档。
    
    此外，我们使用了 [Revise](https://github.com/timholy/Revise.jl) 包来实现代码的热更新，这样在修改代码后不需要重新启动julia kernel，就可以直接在jupyter notebook中看到修改后的效果。

2. 仿真参数配置
    ```julia
    includet(projectdir("config", "example_config.jl"))
    config
    ```
    将参数放置在一个字典中。

3. 定义仿真函数
    ```julia
    function get_E_in_rx_dis_tag_rx(config_now)
   # some code
    end
    ```
4. 仿真
    ```julia
    data, _ = produce_or_load(
            x -> get_E_in_rx_dis_tag_rx(x), config, datadir("example"); verbose = false
        )
    ```
    通过 `produce_or_load` 函数来判断是否需要重新生成数据，如果数据已经存在则直接加载，否则生成新的数据。

注意：Windows的默认文件路径限制约250字符，DrWaston包默认成生成的文件名称可能会超过这一限制导致保存失败，可以用以下方式[解除限制](https://learn.microsoft.com/zh-cn/windows/win32/fileio/maximum-file-path-limitation?tabs=powershell)
如果使用本项目，请引用以下论文：
```bibtex
@article{shang2024liquimager,
  title={LiquImager: Fine-grained Liquid Identification and Container Imaging System with COTS WiFi Devices},
  author={Shang, Fei and Yang, Panlong and Yan, Dawei and Zhang, Sijia and Li, Xiang-Yang},
  journal={Proceedings of the ACM on Interactive, Mobile, Wearable and Ubiquitous Technologies},
  volume={8},
  number={1},
  pages={1--29},
  year={2024},
  publisher={ACM New York, NY, USA}
}

@misc{shang2024limits,
      title={Towards the limits: Sensing Capability Measurement for ISAC Through Channel Encoder}, 
      author={Fei Shang and Haohua Du and Panlong Yang and Xin He and Wen Ma and Xiang-Yang Li},
      year={2024},
      eprint={2405.09497},
      archivePrefix={arXiv},
      primaryClass={cs.IT}
}
```
