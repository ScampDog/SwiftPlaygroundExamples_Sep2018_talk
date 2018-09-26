//: [Previous](@previous)

import Foundation

let d = Double(getJulianDate(month: 4, day: 19, year: 1990) - getJulianDate(month: 12, day: 31, year: 1999))

var omega = 282.9404 + 4.70935e-5*d
var e = 0.016709 - 1.151e-9*d
var M = 356.0470 + 0.9856002585*d
M = adjustAngle(M)

let sun = KeplerElements(type: .elliptical, tp:0.0, aq: 1.0, e: e, i: 0.0, node: 0.0, omega: omega*pi/180.0)

print(sun.description)

//: [Next](@next)
