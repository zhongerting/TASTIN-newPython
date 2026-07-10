# Benke TOPAZ-II single-cell TFE validation data

Generated for the requested validation condition:

- TISA input power: 3412 W row in Benke Appendix B, p.79
- Active-zone power: approximately 3003 W, consistent with P_az = 0.88 P_TISA
- Regulated helium gap pressure: 10 torr
- Source: Benke 1994, Appendix B, Experimental Test Stand Data taken 17-18 August 1994, p.79; also consistent with Table 5-1 example data set.

Important data-quality note:

Benke states that T64 was inoperative and disregarded in the analysis. The Appendix B thermocouple table therefore gives 11 measured collector-sleeve thermocouple values, not 12 measured values. In `benke_sleeve_thermocouple_12pt_digitized.csv`, row thermocouple_index=9 corresponds to T64 and is set to NaN rather than interpolated.

Check values for row TISA=3412 W:

- Mean of all 11 measured T56, T57, T58, T59, T60, T61, T62, T63, T65, T66, T67 = 441.99 °C
- Mean excluding end thermocouples T56 and T67, and excluding inoperative T64 = 453.6789 °C
- Min/Max average = (354.29 + 515.64)/2 = 434.965 °C
- Cooling water inlet = 16.56 °C
- Cooling water outlet #1/#2 = 32.15/32.09 °C
- Average outlet = 32.12 °C = 305.27 K
- Average water temperature rise = 15.56 K
