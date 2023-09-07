# SPHsimulateor

## 基于SPH方法的简单粒子模拟器

本项目中，基于SPH方法实现了2D场景与3D场景下的流体粒子模拟，

### 环境配置

- VisualStudio2022
- glut 3.7?

### 调用方法

src/main.cpp下：

```C++
int main(int argc, char** argv)
{
	SPH_3D(argc, argv,1);
}
```

1. SPH_2D(argc, argv,0)：2D场景下以49个粒子进行模拟
2. SPH_2D(argc, argv,1)：2D场景下以100个粒子进行模拟
3. SPH_3D(argc, argv,0)：3D场景下以27个粒子进行模拟
4. SPH_3D(argc, argv,1)：3D场景下以125个粒子进行模拟

### 运行截图

#### 2D场景下模拟



#### 3D场景下模拟