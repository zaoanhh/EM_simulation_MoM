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