
using PyPlot, Interpolations, DelimitedFiles

function boxcarn(A::AbstractArray, ntosmooth::Int)
           out = zeros(Float64, size(A))
           R = CartesianIndices(A)
           Ifirst, Ilast = first(R), last(R)
           I1 = oneunit(Ifirst)
           for I in R
               n, s = 0, zero(eltype(out))
               for J in max(Ifirst, I-ntosmooth*I1):min(Ilast, I+ntosmooth*I1)
                   s += A[J]
                   n += 1
               end
               out[I] = s/n
           end
           out
end

xy = [-0.0018781493357764578 -1.695087036188732
       0.03190563444800737 -1.2727897388914347
       0.08197434722858449 -0.7486257443884563
       0.13212322491983508 -0.3421037562986724
       0.1756527714154833 -0.22171323866239145
       0.22207970682546951 -0.3532409528172247
       0.2646243701328447 -0.7875343563902888
       0.29499541914796146 -1.3570487860742109
       0.3199267063673842 -1.9437127805771874
       0.3505611543747137 -2.8997652313330287
       0.3673156207054512 -3.4869445716903362
       0.3827187356848373 -4.091015803939533
       0.3940219880897847 -4.678538708199725
       0.3963582226294091 -4.10696289509849
       0.39715987173614276 -1.2833829592304165
       0.3995190105359596 1.2545808520384787
       0.4016720109940449 4.095052679798442
       0.40242785158039396 4.985856619331195
       0.4385364177737059 3.9965357306459
       0.4565964269354099 3.4934722858451672
       0.4774278515803939 2.9233566193311957
       0.5190563444800732 1.8335432890517636
       0.5522789738891433 1.0793346312414114
       0.5854557947778286 0.3923499770957397
       0.624083829592304 -0.2942911131470449
       0.6817796610169491 -0.9629237288135588
       0.7350893266147505 -1.1948579935868064
       0.7800160329821344 -1.1247995877233166
       0.823442510306917 -0.8531550618415018
       0.8654145671094824 -0.44714841960604623
       0.9046037562986715 0.04271644525881957
       0.9424072377462207 0.5661074209803028
       0.97882501145213 1.123024507558406
       1.0030920751259733 1.5111085661933128
       ];

x2 = LinRange(0, 1, 201) # new x-grid
figure()
plot(xy[:,1],xy[:,2])

itp = interpolate((xy[:,1],), xy[:,2], Gridded(Linear()))
yf = boxcarn(itp(x2), 1)
plot(x2, yf)
writedlm("func.txt", yf)
