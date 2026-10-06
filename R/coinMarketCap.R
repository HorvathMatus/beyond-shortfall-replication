# --=========================================================================================
# stiahne zoznam 100 najvacsich cryptocurrencies podla aktualneho marketcapu
getCoinMarketCapList = function() {
  tmp=rjson::fromJSON(file='https://api.coinmarketcap.com/data-api/v3/cryptocurrency/listing?start=1&limit=100&sortBy=market_cap&sortType=desc&cryptoType=all&tagType=all&audited=false')
  if (tmp$status$error_code!=0) stop(paste('Error in parsing coinMarketCap:',tmp$status$error_message))
  data.table::rbindlist(lapply(tmp$data$cryptoCurrencyList,function(x) data.frame(id=x$id,name=x$name,symbol=x$symbol,price=x$quotes[[1]]$price,marketCap=x$quotes[[1]]$marketCap)))
}

# --=========================================================================================
# stiahne dennu historiu danych ids
# ids - vektor idciek, za ktore sa stiahne historia
# timeStart, timeEnd - datumy od-do
# briefOutput - ak true, nebude davat do vystupu timestamp za OHLC samostatne, len jeden den=timestamp
# verbose - ci ma vypisovat progress
# Vylepsenie oproti povodnej verzii: quarter grid sa oreze na [timeStart, timeEnd],
# takze sa nestahuje zbytocne cela historia od 2012.
getCoinMarketCap = function(ids=c(1,1027),timeStart='2012-01-01',timeEnd=Sys.Date()-1,briefOutput=T, verbose=T) {
  timeStart=as.Date(timeStart); timeEnd=as.Date(timeEnd)
  dates=seq.Date(as.Date(cut(timeStart,'quarter')),timeEnd,by='quarter')
  if (dates[length(dates)]<timeEnd) dates=c(dates,timeEnd)
  if (length(dates)<2) dates=c(dates[1],timeEnd)

  df=NULL
  for (id in ids) {
    for (i in 2:length(dates)) {
      tmp=rjson::fromJSON(file=paste0('https://api.coinmarketcap.com/data-api/v3/cryptocurrency/historical?id=',id,'&convertId=2781&timeStart=',as.integer(as.POSIXct(dates[i-1])),'&timeEnd=',as.integer(as.POSIXct(dates[i]))))
      if (length(tmp$data$quotes)>0) {
        if (verbose) cat('downloading from coinMarketCap: ',id,as.character(dates[i]),'\n')
        tmp2=data.table::rbindlist(lapply(tmp$data$quotes, function(x) data.frame(id=id,x)))

        names(tmp2)=gsub('quote.','',names(tmp2))
        tmp2=tmp2[,c('id','timestamp','open','high','low','close','volume','marketCap','timeOpen','timeClose','timeHigh','timeLow')]
        tmp2$timestamp=as.Date(tmp2$timestamp)
        if (briefOutput) tmp2=tmp2[,c('id','timestamp','open','high','low','close','volume','marketCap')]

        df=rbind(df,tmp2)
      }
      Sys.sleep(0.15)  # slusnost voci API
    }
  }

  if (is.null(df)) return(NULL)
  df=df[!duplicated(df),]
  if (sum(duplicated(df[,c('id','timestamp')]))>0) stop('Nieco je zle s coinmarketcap')
  df2=df[(timestamp>=timeStart)&(timestamp<=timeEnd)]
  df2
}
