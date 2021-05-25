# ODE_ellipse
高専時代の卒研内容です。解説はpdfにあります。
Runge-Kutta-Gill methodなどによる惑星シミュレーションをCで書いたやつです。 

#define SOLVER ()の中身を変えると、Euler Methodとかに変更できます。  

(0) Euler  
(1) Heun  
(2) Euler(2nd)   
(3) Runge-Kutta(3rd)  
(4) Runge-Kutta(4th)   
(5) Runge-Kutta-Gill  

e: 離心率  
a: 楕円長半径  

ecact_ellipse: ニュートン力学での惑星運動  
diff_ellipse: Schwarzschils時空における、重力が弱いところでの一般相対性理論の効果を取り入れた時の惑星運動  





