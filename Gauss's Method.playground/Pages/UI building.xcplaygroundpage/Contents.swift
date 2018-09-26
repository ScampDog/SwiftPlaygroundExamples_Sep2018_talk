import UIKit
import PlaygroundSupport

class GaussIteratorViewController: UIViewController {
    var textView : UITextView!
    
    var iterateButton: UIButton!
    var quitButton: UIButton!
    var doneButton: UIButton!
    
    override func loadView() {
        let width = 400
        let height = 300
        let buttonHeight = 50
        let labelHeight = 75
        
        let textHeight = 300
        let view = UIView(frame: CGRect(x: 50, y: 50, width: width, height: height))
        view.backgroundColor = UIColor.black
        
        textView = UITextView(frame: CGRect(x: 20+width/4, y: 5, width: width*3/4, height: textHeight))
        textView.backgroundColor = UIColor.white
        let displayString = "Elliptical Orbit\n\na = 1.4 AU"
        textView.text = displayString
        view.addSubview(textView)
        
        iterateButton = UIButton(frame: CGRect(x: 10, y: 10, width: width/4, height: buttonHeight))
        iterateButton.setTitle("Iterate", for: .normal)
        iterateButton.layer.cornerRadius = 5
        iterateButton.backgroundColor
            = UIColor.green
        iterateButton.tintColor = UIColor.black
        view.addSubview(iterateButton)
        
        doneButton = UIButton(frame: CGRect(x: 10, y: buttonHeight+20, width: width/4, height: buttonHeight))
        doneButton.setTitle("Done", for: .normal)
        doneButton.backgroundColor = UIColor.blue
        view.addSubview(doneButton)
        doneButton.layer.cornerRadius = 5
        
        quitButton = UIButton(frame: CGRect(x: 10, y: 30+2*buttonHeight, width: width/4, height: buttonHeight))
        quitButton.setTitle("Quit", for: .normal)
        quitButton.backgroundColor = UIColor.red
        quitButton.layer.cornerRadius = 5
        view.addSubview(quitButton)
        
        self.view = view
    }
}

let iteratorViewController = GaussIteratorViewController()
PlaygroundPage.current.liveView = iteratorViewController


