# TOA/FOA測位精度評価モデル  

## 概要

このMATLABスクリプトは複数基の衛星によるTOAおよびFOAを複合利用した三次元測位精度を評価するモデルです.  
このモデルでは月の特異重力分布による衛星軌道への影響を再現するために, NASAゴダード宇宙センター(GSFC)が開発したGMAT(General Mission Analalysis Tool)および同NAIF(Navigation and Ancillary Information Facility)が開発したSPICE Toolkitを利用しています.

## 動作確認済み環境

- Windows10 / Windows11
- MATLAB 2024b
- GMAT R2025a

## インストール方法

### MATLAB測位モデルおよび関連ファイルインストール  

``` powershell
git clone https://github.com/TOA_FOA_Navigation_Models.git
```  

### GMAT&SPICE Toolkitインストール

PowerShellで以下のスクリプトを実行する.  

``` powershell
.\tools_installer
```

## 使用方法

### GMATで衛星軌道データを生成する

### MATLABプログラムを実行する  

## ファイル構成

- **constant_files** (シミュレーションに必要な各種天体物理データ群)
  - de421.bsp
  - moon_080317.tf
  - moon_pa_de421_1900-2050.bpc
  - naif0012.tls
  - pck00010.tpc
- **GMAT_Scripts** (GMATシミュレーション用サンプルスクリプト)
  - LEO_PNT_in_gravity.script (重力分布入り低高度円軌道衛星シミュレーション)
  - LEO_PNT_Ephemeris.script (重力分布無し低高度円軌道衛星シミュレーション)
  - ELFO_PNT_in_gravity.script (重力分布入り楕円凍結軌道衛星シミュレーション)
  - ELFO_PNT_Ephemeris.script (重力分布無し楕円凍結軌道衛星シミュレーション)
- toa_and_foa.m   (TOA/FOA逐次測位用スクリプト)
- toa_and_foa_multi_epoch.m   (TOA/FOAマルチエポック測位用スクリプト)
- tools_installer.ps1 (GMAT&SPICE Toolkitインストール用PowerShellスクリプト)
