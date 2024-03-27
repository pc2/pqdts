
function main()
  setprecision(128)
  x=Float64(2.0)
  x2=BigFloat(2.0)
  for i in 1:5000
    x=x^1.01
    x2=x2^1.01
    println(i," ",x," ",x2)
  end
end

main()
