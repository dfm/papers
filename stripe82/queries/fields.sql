SELECT
    p.run,p.camcol,p.field,p.rerun,p.mjd_g,p.ramin,p.ramax,p.decmin,p.decmax
INTO mydb.fields
FROM Stripe82..Field AS p
WHERE p.mjd_g > 0

