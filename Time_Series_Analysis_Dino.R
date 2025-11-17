library(quantmod)
library(zoo)
library(xts)
library(forecast)
library(ggplot2)
library(stats)
library(tseries)
library(seastests)

##############################   PROJEKT  ###################################
################ Modele jednowymiarowych szeregów czasowych #################
###################### na podstawie spółki DINO #############################

start_date <- "2020-01-01"
end_date <- "2024-12-31"

#Pobranie danych finansowych
x <- getSymbols("DNP.WA", src = "yahoo", auto.assign = FALSE, from = start_date, to = end_date)

# Usunięcie braków 
x <- na.omit(x)

# Sprawdzenie, czy dane są szeregiem czasowym
is.ts(x)  
# Wynik jest fałszywy, dlatego trzeba przekształcić dane na szereg czasowy

# Konwersja danych dot. zamknięcia na ts, czyli szereg czasowy
x <- ts(log(x[, "DNP.WA.Close"]), start = c(2020, 1), frequency = 252) 
plot(x, main = "log Dino Stock Prices", ylab = "Price", xlab = "Time")
#wahania rosną wraz z trendem - model jest multiplikatywny

x.ts <- x
frequency(x.ts)
# sezonowość = 252, jest to sezonowość roczna
autoplot(x.ts)

# Wykres, który może ułatwić identyfikację sezonowości
seasonplot(x.ts, col = rainbow(5), main = "Wykres sezonowości", type ='l')
legend("bottomright", legend = 2020:2024, col = rainbow(5), lty = 1)

adf.test(x.ts)
# wynik testu adf wskazuje na to, że szereg jest niestacjonarny, p value > 0.05
acf(x.ts)
#wykres acf to potwierdził

# sprawdzenie stacjonarności po jednym zróżnicowaniu
adf.test(diff(x.ts))
# wynik testu adf wskazuje na to, że szereg jest stacjonarny, p value < 0.05
# szereg jest zintegrowany I(1)

xdiff.ts <- diff(x.ts)

#badanie autokorelacji i autokorelacji cząstkowej
layout(1:2)
acf(xdiff.ts, lag = 504, main = "ACF dla szeregu czasowego")
pacf(xdiff.ts, lag = 504, main = "PACF dla szeregu czasowego")
# pierwsze wnioski:
# Szereg prawdopodobnie nie wykazuje sezonowości.
# Model ARIMA(0,1,1) (lub podobny) może być odpowiedni, 
# biorąc pod uwagę odcięcie ACF i wartości PACF.


# dekompozycja szeregu czasowego    
decomposed_x.ts<-decompose(x.ts,type='multiplicative')
autoplot(decomposed_x.ts)
# Trend rośnie stopniowo do około 2024 roku, po czym wydaje się stabilizować.
# Wzór w sezonowości może sugerować cykliczność związaną z okresem roku lub innymi regularnymi zjawiskami.
# Składnik resztowy nie wskazuje na wyraźne wzorce

###################### Estymacja parametrów ####################################

x_model <- auto.arima(x.ts, trace = TRUE) 
# Najlepszym modelem okazała się SARIMA (0,1,1)(1,0,0)[252]
summary(x_model)

# ręcznie sprawdzimy jeszcze inny model, bez sezonowosci,to może pomóc ocenić
# czy obecność komponentu sezonowego w SARIMA jest uzasadniona, czy też może prowadzić do niepotrzebnej złożoności modelu.
x_model2 <- Arima(x.ts, order = c(0, 1, 1))
summary(x_model2) 

# ręcznie sprawdzimy jeszcze inny model, SARIMA(1,1,1)(0,0,1)[252], który uwzględnia dodatkowy parametr sezonowy MA(1), który może poprawić dopasowanie

x_model3 <- Arima(x.ts, order = c(1,1,1), seasonal = list(order = c(0,0,1), period = 252))
summary(x_model3)

#sprawdzanie reszt
white_noise_test <- Box.test(residuals(x_model), lag = 20, type = "Ljung-Box")
cat("Test białego szumu dla reszt modelu SARIMA (0,1,1)(1,0,0)[252] (Ljung-Box): p-value =", white_noise_test$p.value, "\n")
white_noise_test2 <- Box.test(residuals(x_model2), lag = 20, type = "Ljung-Box")
cat("Test białego szumu dla reszt modelu ARIMA (0,1,1) (Ljung-Box): p-value =", white_noise_test2$p.value, "\n")
white_noise_test2 <- Box.test(residuals(x_model3), lag = 20, type = "Ljung-Box")
cat("Test białego szumu dla reszt modelu SARIMA(1,1,1)(0,0,1)[252] (Ljung-Box): p-value =", white_noise_test2$p.value, "\n")
# we wszystkich przypadkach reszty to biały szum, sugeruje to, że modele dobrze uchwyciły wszystkie struktury w danych

#wykresy dla reszt dla trzech modeli
layout(1:3)
hist(residuals(x_model), main = "Histogram reszt modelu SARIMA (0,1,1)(1,0,0)[252]", xlab = "Reszty")
hist(residuals(x_model2), main = "Histogram reszt modelu ARIMA (0,1,1)", xlab = "Reszty")
hist(residuals(x_model3), main = "Histogram reszt modelu SARIMA(1,1,1)(0,0,1)[252]", xlab = "Reszty")

layout(1:3)
qqnorm(residuals(x_model), main = "Wykres Q-Q reszt modelu SARIMA (0,1,1)(1,0,0)[252]")
qqline(residuals(x_model), col = "red")
qqnorm(residuals(x_model2), main = "Wykres Q-Q reszt modelu ARIMA (0,1,1)")
qqline(residuals(x_model2), col = "red")
qqnorm(residuals(x_model3), main = "Wykres Q-Q reszt modelu SARIMA(1,1,1)(0,0,1)[252]")
qqline(residuals(x_model3), col = "red")

layout(1:3)
acf(residuals(x_model), lag = 252, main = "SARIMA (0,1,1)(1,0,0)[252]")
acf(residuals(x_model2), lag = 252, main = "ARIMA (0,1,1)")
acf(residuals(x_model3), lag = 252, main = "SARIMA(1,1,1)(0,0,1)[252]")

layout(1:3)
pacf(residuals(x_model), lag = 252, main = "SARIMA (0,1,1)(1,0,0)[252]")
pacf(residuals(x_model2), lag = 252, main = "ARIMA (0,1,1)")
pacf(residuals(x_model3), lag = 252, main = "SARIMA(1,1,1)(0,0,1)[252]")

#przedziały ufności dla modelu SARIMA (0,1,1)(1,0,0)[252]

#ma1
-0.0756 + c(-1.96, 1.96) * 0.0284
#sar1
0.0728 + c(-1.96, 1.96) * 0.0346

#przedziały ufności dla modelu ARIMA (0,1,1)

#ma1
-0.0759 + c(-1.96, 1.96) * 0.0284

#przedziały ufności dla trzeciego modelu SARIMA (1,1,1)(0,0,1)[252]

#ar1
0.1481 + c(-1.96, 1.96) * 0.5274
#ma1
-0.2229 + c(-1.96, 1.96) * 0.5224
#sma1
0.0697 + c(-1.96, 1.96) * 0.0337

######################## Prognozowanie #########################################

# prognoza na przyszłość
f <- forecast::forecast(x_model, h = 62)
autoplot(f) 
f2 <- forecast::forecast(x_model2, h = 62)
autoplot(f2) 
f3 <-forecast::forecast(x_model3, h = 62)
autoplot(f3) 


# Wyświetlenie wykresu z danymi, prognozą oraz wykresów ACF i PACF dla modelu SARIMA (0,1,1)(1,0,0)[252]
par(mfrow = c(2, 1))
plot(x.ts, main = "Zamknięcie Dino")
plot(f, main = "Prognoza SARIMA (0,1,1)(1,0,0)[252]", ylab = "Zamknięcie", xlab = "Time")
par(mfrow = c(2,1))
acf_res <- acf(residuals(x_model), lag.max = 20, main = "ACF reszty SARIMA (0,1,1)(1,0,0)[252]")
pacf_res <- pacf(residuals(x_model), lag.max = 20, main = "PACF reszty SARIMA (0,1,1)(1,0,0)[252]")
par(mfrow = c(1, 1))

# Wyświetlenie wykresu z danymi, prognozą oraz wykresów ACF i PACF dla modelu ARIMA (0,1,1)
par(mfrow = c(2, 1))
plot(x.ts, main = "Zamknięcie Dino")
plot(f2, main = "Prognoza ARIMA (0,1,1)", ylab = "Zamknięcie", xlab = "Time")
par(mfrow = c(2, 1))
acf_res <- acf(residuals(x_model2), lag.max = 20, main = "ACF reszty ARIMA (0,1,1)")
pacf_res <- pacf(residuals(x_model2), lag.max = 20, main = "PACF reszty ARIMA (0,1,1)")
par(mfrow = c(1, 1))

# Wyświetlenie wykresu z danymi, prognozą oraz wykresów ACF i PACF dla modelu SARIMA (1,1,1)(0,0,1)[252]
par(mfrow = c(2, 1))
plot(x.ts, main = "Zamknięcie Dino")
plot(f3, main = "Prognoza SARIMA (1,1,1)(0,0,1)[252]", ylab = "Zamknięcie", xlab = "Time")
par(mfrow = c(2, 1))
acf_res <- acf(residuals(x_model3), lag.max = 20, main = "ACF reszty SARIMA (1,1,1)(0,0,1)[252]")
pacf_res <- pacf(residuals(x_model3), lag.max = 20, main = "PACF reszty SARIMA (1,1,1)(0,0,1)[252]")
par(mfrow = c(1, 1))


########## Model uczący dla SARIMA (0,1,1)(1,0,0)[252] ##########################

# Określenie procentowego podziału na zbiór uczący i testowy
procent <- 0.70
p <- procent * round(length(x.ts))
zu <- window(x.ts, end = time(x.ts)[p])
zt <- window(x.ts, start = time(x.ts)[p + 1])

sarima_model_uczacy <- Arima(zu, order = c(0, 1, 1), seasonal = c(1, 0, 0))
summary(sarima_model_uczacy)

#Wygenerowanie prognoz dla zbioru testowego
forecast_values_test <- forecast::forecast(sarima_model_uczacy, h = length(zt))
autoplot(forecast_values_test)

# Test białego szumu dla reszt
white_noise_test <- Box.test(residuals(sarima_model_uczacy), lag = 20, type = "Ljung-Box")
cat("Test białego szumu dla reszt modelu uczącego SARIMA (Ljung-Box): p-value =", white_noise_test$p.value, "\n")

plot(residuals(sarima_model_uczacy), main = "Reszty modelu uczacego SARIMA")
acf_res <- acf(residuals(sarima_model_uczacy), main = "ACF reszt")
pacf_res <- pacf(residuals(sarima_model_uczacy), main = "PACF reszt")


