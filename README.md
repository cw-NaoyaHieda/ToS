# ToS

VaRとESの推定。  
ToS分布を評価分布として用いているが、modeが変化していない問題。

- data csvとRdata
- rmd Rmdファイル置き場(基本、論文の式確認用rmd)
- old_rmd Rmdファイル置き場 ただし、scriptを利用するので、ToSフォルダ直下に移動しないと動かない
- plot plot置き場 原則eps
- script script置き場

ToS/  
├ data/  
│     ├ nikkei225.csv(実験で使用していたデータ)  
│     ├ nky.csv(概論のプロットと投稿論文用に先生に頂いたデータ)  
│     ├ rolling_result.Rdata(ローリング推定の結果 オブジェクトresutlが入っている 特に工夫なし)  
│     └ rolling_result_useolodpara.Rdata(ローリング推定の結果 オブジェクトresutlが入っている 一時点前のパラメータが推定の初期値)  
├ rmd/  
├ old_rmd/
│        ├ 20180209.Rmd(概論のプロットの一部出力用)  
│       └ snow.Rmd(snowパッケージの性能確認rmd)  
├ plot/  
│     ├ ES100_20180209.eps(100~10000での推定ESの挙動)  
│     ├ logreturn_20180209.eps(nky.csvの対数収益率プロット)  
│     ├ MLE_20180209.eps(nky.csvの正規分布とsinh-arcsinh分布の当てはめ)  
│     ├ n225_20180209.eps(nky.csvのプロット)  
│     └ sinh-arcsinh_20180209.eps(sinh-arcsinhのISのプロット)  
├ script/  
│       ├ 20180203.R(概論プロットの0203に出した修正版)  
│       ├ function.R(実験で使用していた、先生が作成した関数)  
│       ├ functions_gaironplot.R(概論の図の作成で使用している関数)  
│       ├ functions_rolling.R(ローリングで使用している関数)  
│       └ gaironplot.R(概論プロットの最終修正版)  
├ rolling.Rmd(ローリング推定を行っているrmd)  
├ toc.css(rmdのcssフォーマット)  
└ ToS.Rproj  


├── README.md  
├── ToS.Rproj  
├── data  
│   ├── nikkei225.csv  
│   ├── nky.csv  
│   ├── rolling_result.Rdata  
│   └── rolling_result_useoldpara.Rdata  
├── old_rmd  
│   ├── 20180209.Rmd  
│   └── snow.Rmd  
├── plot  
│   ├── ES100_20180209.eps  
│   ├── MLE_20180209.eps  
│   ├── logreturn_20180209.eps  
│   ├── n225_20180209.eps  
│   └── sinh-arcsinh_20180209.eps  
├── rmd  
│   ├── A family of skew distributions with mode-invariance through transformation of scale.Rmd  
│   └── Sinh-arcsinh.Rmd  
├── rolling.Rmd  
├── rolling.html  
├── script  
│   ├── 20180203.R  
│   ├── functions.R  
│   ├── functions_gaironplot.R  
│   ├── functions_rolling.R  
│   └── gaironplot.R  
└── toc.css  
