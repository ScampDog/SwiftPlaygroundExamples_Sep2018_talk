import UIKit
import PlaygroundSupport

/// Gauss's method of preliminary orbit determination.
///
/// - Parameters:
///   - tData: Time of observartions
///   - alpha: right ascension data, hours
///   - delta: declination data, degrees
///   - xHigh: high end of range to evaluate Laplace's equation
///   - xLow: low end of range to evaluate Laplace's equation
///   - xGuess: Value for starting Newton's method of solving  Laplace's equation
///   - rCenter: Position of observer at times in tData
///   - mu: dimensionless mass parameter; 1 + mass of seccondary / mass of primary
///   - k: gravitational parameter
///   - ab: factor for abberation corrections, 1/c, au per day
/// - Returns: position and velocity vectors of the body at the middle  of the three observations
// Data for Pallas, page 435-6: times of the observations, right ascension and declination
let tData = [6370.57744, 6378.56789, 6390.65113] // truncated Julian date
let alpha = [  6.38029,   6.40793,   6.38762] // right ascension, hours
let delta = [-24.25104, -26.48060, -29.48400] // declination, degrees
let rCenter = [Vector3D(-0.7735829, -0.5704494, -0.2473703), Vector3D(-0.6780640, -0.6624821, -0.2872733), Vector3D(-0.5091536, -0.7766740, -0.3367798)] // position of Earth at observations

let xHigh = 3.0
let xLow = 2.0
let xGuess = 2.3

var L = Array(repeating: Vector3D(0.0, 0.0, 0.0), count: tData.count)

for iLoop in 0...tData.count-1 {
    let a = alpha[iLoop]/12*pi // convert to radians, from hours
    let d = delta[iLoop]/180*pi // convert to radians, from degrees
    L[iLoop] = Vector3D(cos(d)*cos(a), cos(d)*sin(a), sin(d))
}

// Transform observation times
var tau = [k*(tData[0]-tData[1]), k*(tData[2]-tData[0]), k*(tData[2]-tData[1])]

var D = Array(repeating: [0.0, 0.0, 0.0], count: 3)
for iLoop in 0...2 {
    D[0][iLoop] = dot(rCenter[iLoop], cross(L[1], L[2]))
    D[1][iLoop] = dot(L[0], cross(rCenter[iLoop], L[2]))
    D[2][iLoop] = dot(L[0], cross(L[1], rCenter[iLoop]))
}
let D0 = dot(L[0], cross(L[1], L[2]))
let EE = -2*dot(L[1], rCenter[1])
let FF = dot(rCenter[1], rCenter[1])
let A1 = tau[2]/tau[1]
let B1 = A1*(tau[1]*tau[1] - tau[2]*tau[2])/6
let A3 = -tau[0]/tau[1]
let B3 = A3*(tau[1]*tau[1] - tau[0]*tau[0])/6
let AA = -(A1*D[1][0] - D[1][1] + A3*D[1][2])/D0
//let AA = 1.1292240  // for Rebek-Jewel example
let BB = -(B1*D[1][0] + B3*D[1][2])/D0
//let BB = -1.0828889 // for Rebek-Jewel example

let A = -(AA*AA + AA*EE + FF)
let B = -mu*(2*AA*BB + BB*EE)
let C = -mu*mu*BB*BB

print(String(format: "D0 = %10.7f", D0))
print()
print(String(format: "AA = %10.7f", AA))
print(String(format: "BB = %10.7f", BB))
print()
print(String(format: "EE = %10.7f", EE))
print(String(format: "FF = %10.7f", FF))
print()
print(String(format: " A = %10.7f", A))
print(String(format: " B = %10.7f", B))
print(String(format: " C = %10.7f", C))
print()

let nX = 10
let xSpacing = (xHigh - xLow)/Double(nX)

var x = xLow
for _ in 0...nX {
    let p = AA + mu*BB/(x*x*x)
    let fr = C + x*x*x*(B + x*x*x*(A + x*x))
    print(String(format:"%3.1f   %7.5f   %10.5f", x, p, fr))
    x = x + xSpacing
}
print()

x = xGuess // 2.25 would be a better choice

var dx = 1.0 // bogus value that makes the first check of dx true
while abs(dx) > 1e-7 {
    dx = (C + x*x*x*(B + x*x*x*(A + x*x)))/(x*x*(3*B + x*x*x*(6*A + 8*x*x)))
    x = x - dx
}
print(String(format: "x converges to %10.6f\n", x))

// First guess solution for f and g; next tries will use f and g function
let u2 = mu/(x*x*x)
var F = [1-u2*tau[0]*tau[0]/2, 0, 1-u2*tau[2]*tau[2]/2]
var G = [tau[0]*(1-u2*tau[0]*tau[0]/6), 0.0, tau[2]*(1-u2*tau[2]*tau[2]/6)]

var FG = F[0]*G[2] - F[2]*G[0]
var C_array = [G[2]/FG, -1, -G[0]/FG]
var D_array = [-F[2]/FG, 0, F[0]/FG]

var P = [0.0, 0.0, 0.0]
var P_old = P
var dpMag: Double
var nIterations = 0

var t = [0.0, 0.0, 0.0]
var v: Vector3D

class GaussIterator {
    // Orbital elements: postion, velocity and their time
    var tEpoch = 0.0
    var rElement = Vector3D(0.0, 0.0, 0.0)
    var vElement = Vector3D(0.0, 0.0, 0.0)
    // Observation data
    var L: [Vector3D] // direction, an array of unit vectors
    var t: [Double] // times of the observations
    var rCenter: [Vector3D] // position of observer at those times
    // Variables that need to persist between iterations
    var P: [Double]
    var P_old: [Double]
    var D: [[Double]]
    var D0: Double
    var C_array: [Double]
    var D_array: [Double]
    var dpMag: Double
    var iterationCount = 0
    
    init(L: [Vector3D], t: [Double], rCenter: [Vector3D], C_array: [Double], D_array: [Double], D: [[Double]], D0: Double) {
        self.L = L
        self.t = t
        self.rCenter = rCenter
        P = [0.0, 0.0, 0.0]
        P_old = P
        self.C_array = C_array
        self.D_array = D_array
        self.D = D
        self.D0 = D0
        dpMag = 0.0
    }
    
    func iterate() {
        nIterations += 1
        print("Iteration #\(nIterations)")
        for i in 0...2 {
            let v1 = C_array[0]*D[i][0]
            let v2 = C_array[1]*D[i][1]
            let v3 = C_array[2]*D[i][2]
            P[i] = (v1 + v2 + v3)/(C_array[i]*D0)
            print(String(format:"P(%d)    %9.7f", i, P[i]))
        }
        print()
        var r: [Vector3D]
        r = [P[0]*L[0] - rCenter[0]]
        r.append(P[1]*L[1] - rCenter[1])
        r.append(P[2]*L[2] - rCenter[2])
        
        let v = D_array[0]*r[0] + D_array[2]*r[2]
        
        var dp = [0.0, 0.0, 0.0]
        for i in 0...2 {
            dp[i] = P[i] - P_old[i]
        }
        P_old = P
        t = tData
        for i in 0...2 {
            t[i] = tData[i] - ab*P[i]
        }
        tau[0] = k*(t[0] - t[1])
        tau[2] = k*(t[2] - t[1])
        tau[1] = tau[2] - tau[0]
        
        print(String(format: "t = %12.7f", t[1]))
        print()
        print(String(format: "r[1] = (%10.7f, %10.7f, %10.7f)", r[1].x, r[1].y, r[1].z))
        print(String(format: "v[1] = (%10.7f, %10.7f, %10.7f)", v.x, v.y, v.z))
        print()
        let dpMag = sqrt(dp[0]*dp[0] + dp[1]*dp[1] + dp[2]*dp[2])
        print(String(format: "dp = %10.7f  %10.7f %10.7f; %10.5e", dp[0], dp[1], dp[2], dpMag))
        print()
        let coeffs0 = fgUniversal(r: r[1], v: v, h: tau[0], mu: mu)
        F[0] = (coeffs0.f + F[0])/2
        G[0] = (coeffs0.g + G[0])/2
        let coeffs2 = fgUniversal(r: r[1], v: v, h: tau[2], mu: mu)
        F[2] = (coeffs2.f + F[2])/2
        G[2] = (coeffs2.g + G[2])/2
        FG = F[0]*G[2] - F[2]*G[0]
        C_array = [+G[2]/FG, C_array[1], -G[0]/FG]
        D_array = [-F[2]/FG, D_array[1], +F[0]/FG]
        rElement = r[1]
        vElement = v
    }
}

let iterator = GaussIterator(L: L, t:tData, rCenter:rCenter, C_array: C_array, D_array: D_array, D: D, D0: D0)

class GaussIteratorViewController: UIViewController {
    var gaussIterator: GaussIterator!
    var textView : UITextView!
    
    var iterateButton: UIButton!
    var quitButton: UIButton!
    var doneButton: UIButton!
    
    override func loadView() {
        let width = 400
        let height = 300
        let buttonHeight = 50
        let buttonWidth = width/4
        
        let textHeight = 300
        let view = UIView(frame: CGRect(x: 50, y: 50, width: width, height: height))
        view.backgroundColor = UIColor.black
        
        textView = UITextView(frame: CGRect(x: 20+buttonWidth, y: 10, width: width-buttonWidth+5, height: textHeight))
        textView.backgroundColor = UIColor.white
        let fixedWidthFont = UIFont(name: "Menlo-Regular", size: UIFont.systemFontSize)
        textView.font = fixedWidthFont
        let displayString = "Ready to begin iteration..."
        textView.text = displayString
        view.addSubview(textView)
        
        iterateButton = UIButton(frame: CGRect(x: 10, y: 10, width: buttonWidth, height: buttonHeight))
        iterateButton.setTitle("Iterate", for: .normal)
        iterateButton.layer.cornerRadius = 5
        iterateButton.backgroundColor
            = UIColor(displayP3Red: 0.0, green: 0.75, blue: 0.0, alpha: 1.0)
        iterateButton.tintColor = UIColor.black
        iterateButton.addTarget(self, action: #selector(GaussIteratorViewController.iterate), for: .touchUpInside)
        view.addSubview(iterateButton)
        
        doneButton = UIButton(frame: CGRect(x: 10, y: buttonHeight+20, width: buttonWidth, height: buttonHeight))
        doneButton.setTitle("Accept", for: .normal)
        doneButton.backgroundColor = UIColor.blue
        doneButton.addTarget(self, action: #selector(GaussIteratorViewController.set), for: .touchUpInside)
        view.addSubview(doneButton)
        doneButton.layer.cornerRadius = 5
        
        quitButton = UIButton(frame: CGRect(x: 10, y: 30+2*buttonHeight, width: buttonWidth, height: buttonHeight))
        quitButton.setTitle("Quit", for: .normal)
        quitButton.backgroundColor = UIColor.red
        quitButton.layer.cornerRadius = 5
        quitButton.addTarget(self, action: #selector(GaussIteratorViewController.failed), for: .touchUpInside)
        view.addSubview(quitButton)
        
        self.view = view
    }
    
@objc func iterate() {
    gaussIterator.iterate()
    let r = gaussIterator.rElement
    let v = gaussIterator.vElement
    let t = gaussIterator.tEpoch
    let elements = classicalElements(t0: t, r0: r, v0: v, ec: ec, k: k_earth, mu: 1.0, e_eps: 1e-7)
    textView.text = elements.description
    
    }
    
    @objc func failed() {
        textView.text = "Attempt to initialize orbit failed, try again with new data."
    }
    
    @objc func set() {
        textView.text.append("\n\nSetting elements of orbit with this estimate.")
    }
    
}

let iteratorViewController = GaussIteratorViewController()
iteratorViewController.gaussIterator = iterator
PlaygroundPage.current.liveView = iteratorViewController
