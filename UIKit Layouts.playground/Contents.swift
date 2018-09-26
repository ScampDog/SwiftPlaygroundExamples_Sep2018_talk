import UIKit
import PlaygroundSupport

class ViewController: UIViewController {
    
    var aButton: UIButton!
    override func viewDidLoad() {
        super.viewDidLoad()
        title = "Hello UIKit"
        aButton = UIButton(type: .system)
        aButton.setTitle("My Button", for: .normal)
        aButton.layer.borderColor = UIColor.red.cgColor
        aButton.layer.borderWidth = 2
        aButton.layer.cornerRadius = 5
        view.addSubview(aButton)
        aButton.frame = CGRect(x: 200, y: 200, width: 100, height: 100)
        
        aButton.addTarget(self, action: #selector(ViewController.changeBackgroundColor), for: .touchUpInside)
        
    }
    
    @objc func changeBackgroundColor(){ 
        view.backgroundColor = .blue
    }
}

let controller = ViewController()

PlaygroundPage.current.liveView = UINavigationController(rootViewController: controller)
