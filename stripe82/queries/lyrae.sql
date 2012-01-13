SELECT
    p.objID,p.ra,p.dec,p.g,p.Err_g,p.u,p.r,p.i,p.z
INTO mydb.stars
FROM Stripe82..PhotoPrimary p
WHERE p.ra BETWEEN 0.000000 AND 10.000000
AND p.dec BETWEEN -2.000000 AND 2.000000
AND p.type = 6 AND p.g BETWEEN 18.000000 AND 22.000000
AND p.u - p.g BETWEEN 0.7 AND 1.35
AND p.g - p.r BETWEEN -0.15 AND 0.40
AND p.r - p.i BETWEEN -0.15 AND 0.22
AND p.i - p.z BETWEEN -0.21 AND 0.25

