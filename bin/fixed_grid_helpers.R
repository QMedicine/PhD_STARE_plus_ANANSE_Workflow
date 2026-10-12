fixed_5prime_bins <- function(position, signed_length, width) {
    endpoint0 <- position + ifelse(signed_length < 0L, abs(signed_length)-1L, 0L)-1L
    if (width<=0L || any(!is.finite(endpoint0)) || any(endpoint0<0L)) stop('Invalid 5-prime endpoint or grid width')
    list(endpoint0=endpoint0,bin_start=as.integer(floor(endpoint0/width))*width)
}
